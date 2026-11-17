"""Punjab's burning seasons, 2001-2025, from the NASA FIRMS archive.

Questions:
  1. How did season totals move, and where does the satellite record stop being comparable?
  2. Which districts carry the season, and has that shifted?
  3. Does the season start in the north (Amritsar, Tarn Taran) and finish in the south, every year?

Inputs:  raw_archive/ (fetch_archive.py), GADM 4.1 districts.
Rules:   Oct 1 - Nov 30 only; FIRMS type 0 (vegetation fire) only; VIIRS low-confidence dropped,
         MODIS confidence < 30 dropped; same-day detections within ~400 m (0.004 deg) counted once.
Outputs: archive_summary.json, district_by_year.csv, season_by_year.csv,
         fig_season_totals.png, fig_district_shares.png, fig_burn_timing.png
Run:     python src/season_archive.py
"""
import json, os, pathlib
import numpy as np, pandas as pd, geopandas as gpd
from scipy.stats import kendalltau
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = pathlib.Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
GADM = os.environ.get("GADM_PATH", str(DATA / "gadm41_IND_2.shp"))  # GADM 4.1, India level 2 (see README)
RAW = DATA / "raw_archive"
HERE = ROOT / "outputs"
(HERE / "figures").mkdir(parents=True, exist_ok=True)
(HERE / "tables").mkdir(parents=True, exist_ok=True)
INK, MUTED, LINE, SURFACE = "#0b0b0b", "#52514e", "#8a8880", "#fcfcfb"
FIRE, BLUE, SHADE = "#eb6834", "#2a6fdb", "#f1efe9"
NORTH = ["Amritsar", "Tarn Taran", "Gurdaspur"]
SOUTH = ["Bathinda", "Sangrur", "Mansa", "Barnala", "Muktsar", "Faridkot"]  # GADM spellings checked below

dist = gpd.read_file(GADM)
pb = dist[dist.NAME_1 == "Punjab"][["NAME_2", "geometry"]].to_crs(4326)
names = sorted(pb.NAME_2.unique())
NORTH = [d for d in NORTH if d in names]
SOUTH = [d for d in SOUTH if d in names]


def load(sensor, year):
    f = pd.read_csv(RAW / f"{sensor}_{year}_PBHR.csv")
    f["acq_date"] = pd.to_datetime(f.acq_date)
    f = f[(f.acq_date.dt.month.isin([10, 11])) & (f.type == 0)]
    if sensor == "modis":
        f = f[f.confidence >= 30]
    else:
        f = f[f.confidence != "l"]
    p = gpd.GeoDataFrame(f, geometry=gpd.points_from_xy(f.longitude, f.latitude), crs=4326)
    p = gpd.sjoin(p, pb, predicate="within").drop(columns="index_right")
    key = (p.latitude / 0.004).round().astype(int).astype(str) + "_" + (p.longitude / 0.004).round().astype(int).astype(str) + "_" + p.acq_date.dt.strftime("%Y%m%d")
    p = p.loc[~key.duplicated()].copy()
    p["doy"] = p.acq_date.dt.dayofyear - pd.Timestamp(year=year, month=10, day=1).dayofyear  # 0 = 1 Oct
    p["year"] = year
    return p


viirs = pd.concat([load("viirs-snpp", y) for y in range(2012, 2026)], ignore_index=True)
modis = pd.concat([load("modis", y) for y in range(2001, 2026)], ignore_index=True)
viirs[["latitude", "longitude", "acq_date", "acq_time", "frp", "daynight", "NAME_2", "doy", "year"]].to_csv(HERE / "cache_viirs_punjab_octnov.csv", index=False) if (HERE / "cache").exists() else None

season = pd.DataFrame({
    "viirs_events": viirs.groupby("year").size(),
    "viirs_day_share": viirs.groupby("year").daynight.apply(lambda s: round((s == "D").mean(), 3)),
    "viirs_frp_sum_gw": viirs.groupby("year").frp.sum() / 1000,
}).reindex(range(2001, 2026))
season["modis_events"] = modis.groupby("year").size().reindex(range(2001, 2026))
season.index.name = "year"
season.to_csv(HERE / "tables" / "season_by_year.csv")

dby = viirs.groupby(["year", "NAME_2"]).size().unstack(fill_value=0).reindex(columns=names, fill_value=0)
dby.to_csv(HERE / "tables" / "district_by_year.csv")
share = dby.div(dby.sum(axis=1), axis=0)

# Timing: median day of burning by district group, per year (VIIRS)
timing = viirs.assign(group=np.where(viirs.NAME_2.isin(NORTH), "north", np.where(viirs.NAME_2.isin(SOUTH), "south", "other")))
med = timing[timing.group != "other"].groupby(["year", "group"]).doy.median().unstack()
med["south_minus_north_days"] = med.south - med.north

pre, post = viirs[viirs.year.between(2016, 2021)], viirs[viirs.year.between(2022, 2025)]
tau, p_tau = kendalltau(season.loc[2002:2021, "modis_events"].index, season.loc[2002:2021, "modis_events"].values)
summary = {
    "rules": "Oct-Nov, FIRMS type 0, VIIRS confidence != low, MODIS confidence >= 30, same-day 0.004 deg duplicates removed, inside GADM 4.1 Punjab",
    "viirs_events_by_year": {int(k): int(v) for k, v in season.viirs_events.dropna().items()},
    "modis_events_by_year": {int(k): int(v) for k, v in season.modis_events.dropna().items()},
    "viirs_peak_year": int(season.viirs_events.idxmax()),
    "viirs_2025_vs_2016_2021_mean": round(float(season.loc[2025, "viirs_events"] / season.loc[2016:2021, "viirs_events"].mean()), 3),
    "modis_trend_2002_2021_kendall_tau": round(float(tau), 3), "modis_trend_2002_2021_p": round(float(p_tau), 4),
    "district_share_2016_2021": {k: round(float(v), 3) for k, v in (pre.NAME_2.value_counts(normalize=True)).head(10).items()},
    "district_share_2022_2025": {k: round(float(v), 3) for k, v in (post.NAME_2.value_counts(normalize=True)).head(10).items()},
    "north_districts": NORTH, "south_districts": SOUTH,
    "median_burn_day_from_1_oct": {int(y): {"north": float(r.north), "south": float(r.south), "gap_days": float(r.south_minus_north_days)} for y, r in med.iterrows()},
    "years_south_later_than_north": int((med.south_minus_north_days > 0).sum()), "years_compared": int(med.shape[0]),
    "mean_gap_days": round(float(med.south_minus_north_days.mean()), 1),
}
(HERE / "tables" / "archive_summary.json").write_text(json.dumps(summary, indent=2))

plt.rcParams.update({"font.family": "Helvetica", "axes.edgecolor": LINE, "axes.labelcolor": MUTED, "xtick.color": MUTED, "ytick.color": MUTED})

# Fig 1: season totals
fig, ax = plt.subplots(figsize=(10, 5.2), dpi=200)
ax.axvspan(2021.5, 2025.5, color=SHADE, zorder=0)
ax.text(2023.5, 0.97, "From 2022 many fires are lit\nafter the 1:30 pm satellite pass,\nso counts here undercount burning", transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=8.5, color=MUTED)
ax.plot(season.index, season.viirs_events, color=FIRE, lw=2.2, marker="o", ms=4, label="VIIRS 375 m (2012-2025)")
ax.plot(season.index, season.modis_events, color=BLUE, lw=2.2, marker="o", ms=4, label="MODIS 1 km (2001-2025)")
ax.set_ylabel("Fire detections, 1 Oct - 30 Nov, Punjab")
ax.set_xlim(2000.5, 2025.5)
ax.spines[["top", "right"]].set_visible(False)
ax.legend(frameon=False, loc="upper left", fontsize=9)
ax.set_title("Twenty-five burning seasons in Punjab, as the satellites saw them", loc="left", fontsize=13, color=INK, fontweight="bold")
fig.text(0.01, 0.01, "Source: NASA FIRMS archive (VIIRS S-NPP, MODIS C6.1), vegetation fires only, duplicates removed. Shift: Jethva 2026; Misra et al. 2026.", fontsize=7, color=MUTED)
fig.tight_layout(rect=(0, 0.03, 1, 1))
fig.savefig(HERE / "figures" / "fig_season_totals.png")

# Fig 2: district shares heatmap (top 12 districts by total)
top = dby.sum().sort_values(ascending=False).head(12).index
fig, ax = plt.subplots(figsize=(10, 5.6), dpi=200)
m = share[top].T.values
im = ax.imshow(m, aspect="auto", cmap="Oranges", vmin=0, vmax=max(0.2, float(m.max())))
ax.set_yticks(range(len(top)), top, fontsize=9)
ax.set_xticks(range(len(share.index)), share.index, fontsize=8, rotation=0)
for i in range(m.shape[0]):
    for j in range(m.shape[1]):
        ax.text(j, i, f"{m[i, j]*100:.0f}", ha="center", va="center", fontsize=6.5, color=INK if m[i, j] < 0.12 else "white")
ax.set_title("Each district's share of Punjab's season (%, VIIRS)", loc="left", fontsize=13, color=INK, fontweight="bold")
for s in ax.spines.values():
    s.set_visible(False)
fig.text(0.01, 0.01, "Source: NASA FIRMS VIIRS S-NPP archive, 1 Oct - 30 Nov each year; GADM 4.1 district lines (Malerkotla, created 2021, sits inside Sangrur here).", fontsize=7, color=MUTED)
fig.tight_layout(rect=(0, 0.03, 1, 1))
fig.savefig(HERE / "figures" / "fig_district_shares.png")

# Fig 3: timing, north vs south
fig, ax = plt.subplots(figsize=(10, 4.8), dpi=200)
ax.plot(med.index, med.north, color=BLUE, lw=2.2, marker="o", ms=4, label="North: " + ", ".join(NORTH))
ax.plot(med.index, med.south, color=FIRE, lw=2.2, marker="o", ms=4, label="South: " + ", ".join(SOUTH))
ax.set_yticks([10, 20, 30, 40, 50], ["11 Oct", "21 Oct", "31 Oct", "10 Nov", "20 Nov"])
ax.set_ylabel("Median day of burning")
ax.spines[["top", "right"]].set_visible(False)
ax.legend(frameon=False, loc="upper left", fontsize=8.5)
ax.set_title(f"The south burns later than the north in {summary['years_south_later_than_north']} of {summary['years_compared']} seasons", loc="left", fontsize=13, color=INK, fontweight="bold")
fig.text(0.01, 0.01, "Source: NASA FIRMS VIIRS S-NPP archive, 2012-2025; median detection date within each district group.", fontsize=7, color=MUTED)
fig.tight_layout(rect=(0, 0.03, 1, 1))
fig.savefig(HERE / "figures" / "fig_burn_timing.png")
print(json.dumps({k: summary[k] for k in ["viirs_peak_year", "viirs_2025_vs_2016_2021_mean", "modis_trend_2002_2021_kendall_tau", "modis_trend_2002_2021_p", "years_south_later_than_north", "years_compared", "mean_gap_days"]}, indent=1))
print(season.to_string())
print(med.round(1).to_string())
print("pre", summary["district_share_2016_2021"]); print("post", summary["district_share_2022_2025"])
