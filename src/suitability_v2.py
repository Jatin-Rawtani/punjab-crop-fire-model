"""Suitability model v2: where should straw collection go in Punjab?

Rebuilds my May 2026 ArcGIS model (an ArcGIS StoryMap) in open-source Python, for all of Punjab, with:
  - six seasons of fire data (2016-2021, VIIRS S-NPP) instead of 40 points; 2022-2025 held out for checking
  - a 10 m land-cover crop mask (Impact Observatory / Esri LULC 2023, via Microsoft Planetary Computer)
  - real straw buyers (PEDA's 11 paddy-straw power plants + 6 operating bio-CNG plants), coal plants as a variant
  - OpenStreetMap roads (motorway to tertiary), distance in km
  - straw supply inside a 15 km collection radius
  - 2,000 Monte Carlo weight draws (Ligmann-Zielinska et al. 2014) so zones are reported with how robust they are

Two questions, two scenarios:
  A. Feed existing buyers: collection centres near current buyers.
  B. Site a new buyer: where burning is high and no buyer already draws straw (gap areas).

Inputs:  data/raw_archive/ (fetch_archive.py), data/buyers.geojson, data/northern-zone-latest.osm.pbf, GADM 4.1
Outputs: outputs/model/ (rasters, zones_A.geojson, zones_B.geojson, summary_v2.json, figures, web overlays)
Run:     python src/suitability_v2.py
"""
import json, os, pathlib
import numpy as np, pandas as pd, geopandas as gpd, rasterio
from rasterio.features import rasterize
from rasterio.transform import from_origin
from rasterio.warp import reproject, Resampling, calculate_default_transform
from scipy import ndimage
from scipy.stats import spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patheffects as pe

ROOT = pathlib.Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
OUT = ROOT / "outputs" / "model"
OUT.mkdir(parents=True, exist_ok=True)
GADM = os.environ.get("GADM_PATH", str(DATA / "gadm41_IND_2.shp"))  # GADM 4.1, India level 2 (see README)
PBF = DATA / "northern-zone-latest.osm.pbf"  # Geofabrik, India northern zone (see README)
CRS = "EPSG:32643"
RES = 1000.0
RNG = np.random.default_rng(42)
N_DRAWS = 2000
RADIUS_KM = 15            # collection radius (Kurinji & Kumar 2021, as cited in the May StoryMap; verify)
BUYER_FULL_KM, BUYER_ZERO_KM = 50, 100   # within 50 km = fully viable; beyond 100 km = not viable (assumption)
STRAW_T_PER_HA = 4.81     # baled straw, ex-situ adopters: 19.48 q/acre (Bhullar et al. 2024)
PADDY_SHARE = 0.75        # share of kharif cropland under paddy (assumption; replace with district statistics)

# ---------- grid ----------
dist = gpd.read_file(GADM)
pb = dist[dist.NAME_1 == "Punjab"][["NAME_2", "geometry"]].to_crs(CRS)
state = pb.dissolve().geometry.iloc[0]
minx, miny, maxx, maxy = [np.floor(v / RES) * RES for v in state.bounds[:2]] + [np.ceil(v / RES) * RES for v in state.bounds[2:]]
W, H = int((maxx - minx) / RES), int((maxy - miny) / RES)
T = from_origin(minx, maxy, RES, RES)
inside = rasterize([(state, 1)], out_shape=(H, W), transform=T, fill=0).astype(bool)
district_id = rasterize([(g, i + 1) for i, g in enumerate(pb.geometry)], out_shape=(H, W), transform=T, fill=0)
print("grid", W, H, "cells in Punjab", int(inside.sum()))


def save(name, arr, dtype="float32"):
    with rasterio.open(OUT / f"{name}.tif", "w", driver="GTiff", height=H, width=W, count=1, dtype=dtype, crs=CRS, transform=T, nodata=-9999, compress="deflate") as d:
        a = np.where(inside, arr, -9999).astype(dtype)
        d.write(a, 1)


# ---------- land cover (Impact Observatory 10 m LULC 2023, class 5 = crops, 7 = built) ----------
lc_cache = OUT / "landcover_frac.npz"
if lc_cache.exists():
    z = np.load(lc_cache); crop, built = z["crop"], z["built"]
else:
    import pystac_client, planetary_computer, odc.stac
    cat = pystac_client.Client.open("https://planetarycomputer.microsoft.com/api/stac/v1", modifier=planetary_computer.sign_inplace)
    bbox4326 = tuple(pb.to_crs(4326).total_bounds)
    items = cat.search(collections=["io-lulc-annual-v02"], bbox=bbox4326, datetime="2023").item_collection()
    print("LULC tiles", len(items))
    ds = odc.stac.load(items, bands=["data"], crs=CRS, resolution=100, bbox=bbox4326, resampling="nearest", groupby="solar_day")
    lc = ds["data"].isel(time=0).values
    src_T = ds.odc.geobox.affine
    crop, built = np.zeros((H, W), "float32"), np.zeros((H, W), "float32")
    for cls, dst in [(5, crop), (7, built)]:
        reproject((lc == cls).astype("float32"), dst, src_transform=src_T, src_crs=CRS, dst_transform=T, dst_crs=CRS, resampling=Resampling.average)
    np.savez_compressed(lc_cache, crop=crop, built=built)
eligible = inside & (crop >= 0.5) & (built < 0.3)
print("eligible crop cells", int(eligible.sum()), "mean crop frac in Punjab", round(float(crop[inside].mean()), 3))

# ---------- fires ----------
def fires(years):
    frames = []
    for y in years:
        f = pd.read_csv(DATA / "raw_archive" / f"viirs-snpp_{y}_PBHR.csv")
        f["acq_date"] = pd.to_datetime(f.acq_date)
        f = f[f.acq_date.dt.month.isin([10, 11]) & (f.type == 0) & (f.confidence != "l")]
        key = (f.latitude / 0.004).round().astype(int).astype(str) + "_" + (f.longitude / 0.004).round().astype(int).astype(str) + "_" + f.acq_date.dt.strftime("%Y%m%d")
        frames.append(f.loc[~key.duplicated()])
    f = pd.concat(frames)
    p = gpd.GeoDataFrame(f, geometry=gpd.points_from_xy(f.longitude, f.latitude), crs=4326).to_crs(CRS)
    col = ((p.geometry.x - minx) // RES).astype(int); row = ((maxy - p.geometry.y) // RES).astype(int)
    ok = (col >= 0) & (col < W) & (row >= 0) & (row < H)
    grid = np.zeros((H, W), "float32")
    np.add.at(grid, (row[ok].values, col[ok].values), 1)
    return grid / len(years), p


pre_counts, pre_pts = fires(range(2016, 2022))
post_counts, post_pts = fires(range(2022, 2026))
burn_pre = ndimage.gaussian_filter(pre_counts, sigma=2)     # detections / km2 / season, smoothed over ~2 km
burn_post = ndimage.gaussian_filter(post_counts, sigma=2)
rho, _ = spearmanr(burn_pre[eligible], burn_post[eligible])

# ---------- supply within 15 km ----------
r = RADIUS_KM
yy, xx = np.mgrid[-r:r + 1, -r:r + 1]
disk = ((xx ** 2 + yy ** 2) <= r ** 2).astype("float32")
crop_km2 = np.where(inside, crop, 0)
supply_km2 = ndimage.convolve(crop_km2, disk, mode="constant")   # km2 of cropland within 15 km

# ---------- roads ----------
road_cache = OUT / "roads_punjab.gpkg"
if road_cache.exists():
    roads = gpd.read_file(road_cache)
else:
    import pyogrio
    b = tuple(pb.to_crs(4326).total_bounds)
    roads = pyogrio.read_dataframe(PBF, layer="lines", bbox=b, where="highway IN ('motorway','trunk','primary','secondary','tertiary','motorway_link','trunk_link','primary_link','secondary_link')")
    roads = roads[["highway", "name", "geometry"]].to_crs(CRS)
    roads = roads[roads.intersects(state.buffer(5000))]
    roads.to_file(road_cache, driver="GPKG")
road_r = rasterize([(g, 1) for g in roads.geometry], out_shape=(H, W), transform=T, fill=0, all_touched=True)
road_km = ndimage.distance_transform_edt(road_r == 0) * RES / 1000

# ---------- buyers ----------
buyers = gpd.read_file(DATA / "buyers.geojson").to_crs(CRS)


def dist_km(gdf):
    m = np.zeros((H, W), "uint8")
    col = ((gdf.geometry.x - minx) // RES).astype(int).clip(0, W - 1); row = ((maxy - gdf.geometry.y) // RES).astype(int).clip(0, H - 1)
    m[row.values, col.values] = 1
    d, idx = ndimage.distance_transform_edt(m == 0, return_indices=True)
    return d * RES / 1000, idx


straw_buyers = buyers[buyers["type"].isin(["biomass_power", "cbg"])]
buyer_km, buyer_idx = dist_km(straw_buyers)
buyer_km_coal, _ = dist_km(buyers)

# ---------- memberships (0-1) ----------
def pct_rank(a):
    v = a[eligible]
    ranks = np.empty_like(v); ranks[np.argsort(v, kind="stable")] = np.arange(v.size)
    out = np.zeros_like(a, dtype="float32"); out[eligible] = ranks / (v.size - 1)
    return out


def ramp_down(d, full, zero):
    return np.clip((zero - d) / (zero - full), 0, 1).astype("float32")


f_burn = pct_rank(burn_pre)
f_supply = pct_rank(supply_km2)
f_road = ramp_down(road_km, 2, 10)                       # <= 2 km full, >= 10 km zero (assumption)
f_market = ramp_down(buyer_km, BUYER_FULL_KM, BUYER_ZERO_KM)
f_market_coal = ramp_down(buyer_km_coal, BUYER_FULL_KM, BUYER_ZERO_KM)
f_gap = np.clip((buyer_km - RADIUS_KM) / (50 - RADIUS_KM), 0, 1).astype("float32")  # 0 inside a buyer's 15 km pull, 1 beyond 50 km

for n, a in [("burn_pre", burn_pre), ("burn_post", burn_post), ("crop_frac", crop), ("supply_km2_15km", supply_km2), ("road_km", road_km), ("buyer_km", buyer_km)]:
    save(n, a)

BASE = np.array([0.40, 0.30, 0.20, 0.10])  # burning, market (or gap), road, supply: the May 2026 weights, cropland -> supply


def run(name, layers):
    stack = np.stack([l[eligible] for l in layers])               # 4 x n
    base_score = BASE @ stack
    draws = RNG.dirichlet(BASE * 20, size=N_DRAWS)                # mean = base weights, roughly +-0.1
    top = np.zeros(stack.shape[1], "float32"); acc = np.zeros(stack.shape[1], "float64"); acc2 = np.zeros_like(acc)
    for w in draws:
        s = w @ stack
        top += s >= np.quantile(s, 0.9)
        acc += s; acc2 += s ** 2
    p_top = top / N_DRAWS
    sd = np.sqrt(acc2 / N_DRAWS - (acc / N_DRAWS) ** 2)
    full = lambda v: np.where(eligible, 0, np.nan).astype("float32") if v is None else _put(v)
    score_r, ptop_r, sd_r = _put(base_score), _put(p_top), _put(sd)
    save(f"score_{name}", score_r); save(f"ptop10_{name}", ptop_r); save(f"sd_{name}", sd_r)
    robust = np.nan_to_num(ptop_r) >= 0.8
    # Hubs: greedy pick of the highest-scoring robust cell, then block 2 x radius so catchments do not overlap.
    pre_15 = ndimage.convolve(np.where(inside, pre_counts, 0), disk, mode="constant")    # detections / season within 15 km
    post_15 = ndimage.convolve(np.where(inside, post_counts, 0), disk, mode="constant")
    cand = np.where(robust, np.nan_to_num(score_r, nan=-1), -1).copy()
    rows_, cols_ = np.mgrid[0:H, 0:W]
    zones = []
    while cand.max() > 0:
        rr, cc = np.unravel_index(np.argmax(cand), cand.shape)
        near = (rows_ - rr) ** 2 + (cols_ - cc) ** 2 <= (2 * RADIUS_KM) ** 2
        catch = (rows_ - rr) ** 2 + (cols_ - cc) ** 2 <= RADIUS_KM ** 2
        robust_in = int((robust & catch).sum())
        cand[near] = -1
        if robust_in < 25:            # a hub needs a real block of robust land around it, not a lone cell
            continue
        x, y = minx + (cc + 0.5) * RES, maxy - (rr + 0.5) * RES
        did = district_id[rr, cc]
        bi = buyer_idx[:, rr, cc]
        bx, by = minx + (bi[1] + 0.5) * RES, maxy - (bi[0] + 0.5) * RES
        nb = straw_buyers.iloc[((straw_buyers.geometry.x - bx) ** 2 + (straw_buyers.geometry.y - by) ** 2).values.argmin()]
        sup = float(supply_km2[rr, cc])
        zones.append(dict(zone=len(zones) + 1, district=pb.NAME_2.iloc[did - 1] if did else "", score=round(float(score_r[rr, cc]), 3),
                          p_top10=round(float(ptop_r[rr, cc]), 2), robust_km2_in_catchment=robust_in,
                          fires_per_season_15km_2016_21=int(pre_15[rr, cc]), fires_per_season_15km_2022_25=int(post_15[rr, cc]),
                          crop_km2_within_15km=round(sup), straw_t_within_15km=int(sup * 100 * PADDY_SHARE * STRAW_T_PER_HA),
                          nearest_buyer=nb["name"], nearest_buyer_km=round(float(buyer_km[rr, cc]), 1), road_km=round(float(road_km[rr, cc]), 1), x=x, y=y))
    z = pd.DataFrame(zones).sort_values(["fires_per_season_15km_2016_21"], ascending=False).reset_index(drop=True)
    z.insert(0, "rank", range(1, len(z) + 1))
    gz = gpd.GeoDataFrame(z, geometry=gpd.points_from_xy(z.x, z.y), crs=CRS)
    gz.to_crs(4326).to_file(OUT / f"zones_{name}.geojson", driver="GeoJSON")
    z.drop(columns=["x", "y"]).to_csv(OUT / f"zones_{name}.csv", index=False)
    top_cells = _put((base_score >= np.quantile(base_score, 0.9)).astype("float32"))
    lift = float(np.nansum(np.where(top_cells == 1, post_counts, 0)) / np.nansum(np.where(eligible, post_counts, 0)) / 0.10)
    return dict(robust_km2=int(robust.sum()), hubs=int(len(z)), lift_post_2022=round(lift, 2), top_zones=z.head(10).drop(columns=["x", "y"]).to_dict("records")), ptop_r


def _put(v):
    out = np.full((H, W), np.nan, "float32"); out[eligible] = v; return out


resA, pA = run("A", [f_burn, f_market, f_road, f_supply])
resA_coal, _ = run("A_with_coal", [f_burn, f_market_coal, f_road, f_supply])
resB, pB = run("B", [f_burn, f_gap, f_road, f_supply])

summary = {
    "grid": f"{W} x {H} at 1 km, EPSG:32643; {int(inside.sum())} km2 Punjab, {int(eligible.sum())} km2 eligible (crop >= 50%, built < 30%)",
    "fires_2016_2021_detections": int(len(pre_pts)), "fires_2022_2025_detections": int(len(post_pts)),
    "spatial_persistence_spearman_pre_vs_post": round(float(rho), 3),
    "share_eligible_within_50km_of_straw_buyer": round(float((buyer_km[eligible] <= 50).mean()), 3),
    "share_eligible_within_15km_of_straw_buyer": round(float((buyer_km[eligible] <= 15).mean()), 3),
    "straw_buyer_demand_t_yr_power_plants": int(straw_buyers.straw_t_yr.fillna(0).sum()),
    "assumptions": {"radius_km": RADIUS_KM, "buyer_full_km": BUYER_FULL_KM, "buyer_zero_km": BUYER_ZERO_KM, "road_full_km": 2, "road_zero_km": 10,
                    "straw_t_per_ha": STRAW_T_PER_HA, "paddy_share": PADDY_SHARE, "weights_base": BASE.tolist(), "dirichlet_concentration": 20, "draws": N_DRAWS, "robust_rule": "P(top 10%) >= 0.8"},
    "scenario_A_feed_existing_buyers": resA, "scenario_A_with_coal": {k: v for k, v in resA_coal.items() if k != "top_zones"}, "scenario_B_site_new_buyer": resB,
}
(OUT / "summary_v2.json").write_text(json.dumps(summary, indent=2, default=str))

# ---------- figures ----------
INK, MUTED = "#0b0b0b", "#52514e"
cmap = LinearSegmentedColormap.from_list("rob", ["#f4f1ea", "#f9c58d", "#eb6834", "#8f2d0b"])
ext = (minx, maxx, miny, maxy)
for name, P, res, title in [("A", pA, resA, "Collection centres to feed today's buyers"), ("B", pB, resB, "Where a new straw buyer would face the least competition")]:
    fig, ax = plt.subplots(figsize=(9, 9.6), dpi=200)
    base = np.where(inside, 0.0, np.nan)
    ax.imshow(np.where(inside & ~eligible, 1, np.nan), extent=ext, cmap=LinearSegmentedColormap.from_list("g", ["#e7e5df", "#e7e5df"]), zorder=1)
    ax.imshow(P, extent=ext, cmap=cmap, vmin=0, vmax=1, zorder=2)
    pb.boundary.plot(ax=ax, color="#8a8880", lw=0.4, zorder=3)
    for kind, mk, col, lab in [("biomass_power", "s", "#2a6fdb", "Paddy-straw power plant (PEDA)"), ("cbg", "^", "#1b8a5a", "Bio-CNG plant (approx. for 5 of 6)"), ("coal_cofiring", "x", INK, "Coal plant (pellet co-firing)")]:
        b = buyers[buyers["type"] == kind]
        ax.scatter(b.geometry.x, b.geometry.y, marker=mk, s=38, c=col, edgecolors="white", linewidths=0.6, zorder=5, label=lab)
    z = gpd.read_file(OUT / f"zones_{name}.geojson").to_crs(CRS).head(5)
    for _, rr in z.iterrows():
        ax.annotate(f"{rr['rank']}", (rr.geometry.x, rr.geometry.y), fontsize=9, fontweight="bold", color="white", ha="center", va="center", zorder=6,
                    bbox=dict(boxstyle="circle,pad=0.25", fc=INK, ec="none"))
    for dn in ["Amritsar", "Tarn Taran", "Firozpur", "Fazilka", "Muktsar", "Bathinda", "Mansa", "Sangrur", "Barnala", "Moga", "Ludhiana", "Patiala", "Jalandhar", "Hoshiarpur", "Gurdaspur", "Faridkot"]:
        g = pb[pb.NAME_2 == dn]
        if len(g):
            rp = g.geometry.iloc[0].representative_point()
            ax.text(rp.x, rp.y + 9000, dn, fontsize=6.5, color="#52514e", ha="center", zorder=4, style="italic", path_effects=[pe.withStroke(linewidth=2, foreground="white")])
    ax.set_axis_off()
    ax.legend(loc="lower left", frameon=False, fontsize=8)
    ax.set_title(title, loc="left", fontsize=13, fontweight="bold", color=INK)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, 1))
    cb = fig.colorbar(sm, ax=ax, shrink=0.35, pad=0.01, location="right")
    cb.set_label("Share of 2,000 weightings that put this km² in the top 10%", fontsize=8, color=MUTED)
    cb.outline.set_visible(False)
    fig.text(0.02, 0.015, "Fires: NASA FIRMS VIIRS 2016-2021 · Crops: Impact Observatory/Esri LULC 2023 · Roads: OpenStreetMap · Buyers: PEDA, WRI GPPD · Districts: GADM 4.1\n"
             "Numbered: candidate hubs (15 km catchments, no overlap), ranked by burning inside the catchment. Grey: not cropland.", fontsize=6.5, color=MUTED)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(OUT / f"fig_suitability_{name}.png")
    plt.close(fig)

# ---------- web overlays (EPSG:3857 PNGs for the story map) ----------
import matplotlib.image as mpimg
for name, P in [("A", pA), ("B", pB), ("burn_pre", np.where(eligible, f_burn, np.nan))]:
    dst_crs = "EPSG:3857"
    t2, w2, h2 = calculate_default_transform(CRS, dst_crs, W, H, minx, miny, maxx, maxy)
    dst = np.full((h2, w2), np.nan, "float32")
    reproject(P.astype("float32"), dst, src_transform=T, src_crs=CRS, dst_transform=t2, dst_crs=dst_crs, resampling=Resampling.nearest, src_nodata=np.nan, dst_nodata=np.nan)
    v = np.nan_to_num(dst, nan=0); rgba = cmap(v); rgba[..., 3] = np.where(np.isnan(dst) | (v < 0.03), 0, np.clip(0.25 + 0.7 * v, 0, 0.92))
    mpimg.imsave(OUT / f"web_{name}.png", rgba)
    l, b2, r2, t3 = t2.c, t2.f + t2.e * h2, t2.c + t2.a * w2, t2.f
    from pyproj import Transformer
    tr = Transformer.from_crs(dst_crs, 4326, always_xy=True)
    (lon0, lat0), (lon1, lat1) = tr.transform(l, b2), tr.transform(r2, t3)
    json.dump({"south": lat0, "west": lon0, "north": lat1, "east": lon1}, open(OUT / f"web_{name}_bounds.json", "w"))

print(json.dumps({k: v for k, v in summary.items() if not k.startswith("scenario")}, indent=1))
for k in ["scenario_A_feed_existing_buyers", "scenario_A_with_coal", "scenario_B_site_new_buyer"]:
    s = summary[k]; print(k, {kk: vv for kk, vv in s.items() if kk != "top_zones"})
    for zz in s.get("top_zones", [])[:10]:
        print("   ", zz["rank"], zz["district"], "robust km2", zz["robust_km2_in_catchment"], "p", zz["p_top10"], "fires/season", zz["fires_per_season_15km_2016_21"], "->", zz["fires_per_season_15km_2022_25"], "buyer", zz["nearest_buyer"], zz["nearest_buyer_km"], "road", zz["road_km"])
