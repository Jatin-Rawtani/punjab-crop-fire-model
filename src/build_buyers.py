"""Straw buyer layer for Punjab (who pays for straw, and where).

Sources (read 8 Nov 2026):
  PEDA, Biomass Power Projects, Table I (11 plants, 101.5 MW, 11.7 lakh t straw/yr): https://www.peda.gov.in/biomass-power-projects.php
  PEDA, Bio-CNG/CBG Projects (57 allotted, 6 operational = 107.48 TPD): https://www.peda.gov.in/bio-cbg-compressed-piogas-projects.php
  Operational CBG districts (Khanna, Sangrur, Moga, Fatehgarh Sahib, Patiala, Hoshiarpur): BioEnergy Times / news reports; only
     Bhuttal Kalan (Sangrur) has a village-level location, the other five sit at the town/district centre (flagged 'approx').
  Coal plants (5% biomass co-firing rule; buy pellets, not loose straw): WRI Global Power Plant Database + OpenStreetMap.
Geocoding: OpenStreetMap Nominatim, village + tehsil + district, cached in geocode_cache.json. Topic-only queries.
Output: buyers.geojson, buyers.csv
"""
import json, pathlib, time, urllib.parse, urllib.request
import geopandas as gpd, pandas as pd

HERE = pathlib.Path(__file__).resolve().parent.parent / "data"
CACHE = HERE / "geocode_cache.json"
cache = json.loads(CACHE.read_text()) if CACHE.exists() else {}

ROWS = [  # name, type, place query, capacity, straw t/yr, precision
    ("Malwa Power (6 MW)", "biomass_power", "Gulabewala, Sri Muktsar Sahib, Punjab", 6.0, 92000),
    ("Dee Development Engineers (8 MW)", "biomass_power", "Gadda Dhob, Abohar, Fazilka, Punjab", 8.0, 95000),
    ("Universal Biomass Energy (14.5 MW)", "biomass_power", "Channu, Malout, Sri Muktsar Sahib, Punjab", 14.5, 170000),
    ("Green Planet Energy, Binjon (6 + 4 MW)", "biomass_power", "Binjon, Garhshankar, Hoshiarpur, Punjab", 10.0, 85000),
    ("Green Planet Energy, Bir Pind (6 MW)", "biomass_power", "Bir Pind, Nakodar, Jalandhar, Punjab", 6.0, 72000),
    ("Viaton Energy (10 MW)", "biomass_power", "Khokhar Khurd, Mansa, Punjab", 10.0, 110000),
    ("Sampuran Agri Venture (1 MW)", "biomass_power", "Panchewali, Fazilka, Punjab", 1.0, 6000),
    ("Sukhbir Agro Energy, Jaitu (18 MW)", "biomass_power", "Sedda Singh Wala, Jaitu, Faridkot, Punjab", 18.0, 220000),
    ("Sukhbir Agro Energy, Ferozepur (18 MW)", "biomass_power", "Hukumat Singh Wala, Firozpur, Punjab", 18.0, 220000),
    ("Sukhbir Agro Energy, Jalkheri (10 MW)", "biomass_power", "Jalkheri, Fatehgarh Sahib, Punjab", 10.0, 100000),
    ("CBG plant, Bhuttal Kalan (33.23 TPD)", "cbg", "Bhutal Kalan, Lehragaga, Sangrur, Punjab", 33.23, None),
    ("CBG plant, Khanna (approx)", "cbg", "Khanna, Ludhiana, Punjab", None, None),
    ("CBG plant, Moga district (approx)", "cbg", "Moga, Punjab", None, None),
    ("CBG plant, Fatehgarh Sahib district (approx)", "cbg", "Fatehgarh Sahib, Punjab", None, None),
    ("CBG plant, Patiala district (approx)", "cbg", "Patiala, Punjab", None, None),
    ("CBG plant, Hoshiarpur district (approx)", "cbg", "Hoshiarpur, Punjab", None, None),
]
COAL = [  # from WRI GPPD (lat, lon), Punjab only, operating
    ("GGSSTP Ropar (840 MW)", 31.0417, 76.5847), ("GHTP Lehra Mohabbat (920 MW)", 30.2666, 75.1661),
    ("Goindwal Sahib TPP (540 MW)", 31.3841, 75.1464), ("Rajpura TPP (1,400 MW)", 30.5570, 76.5770),
    ("Talwandi Sabo TPP (1,980 MW)", 29.9240, 75.2372),
]


def geocode(q):
    if q in cache:
        return cache[q]
    url = "https://nominatim.openstreetmap.org/search?" + urllib.parse.urlencode({"q": q, "format": "json", "limit": 1, "countrycodes": "in"})
    req = urllib.request.Request(url, headers={"User-Agent": "punjab-crop-fire-model/1.0"})
    res = json.loads(urllib.request.urlopen(req, timeout=30).read())
    time.sleep(1.1)
    cache[q] = [float(res[0]["lat"]), float(res[0]["lon"]), res[0]["display_name"]] if res else None
    CACHE.write_text(json.dumps(cache, indent=1))
    return cache[q]


out = []
for name, kind, q, cap, straw in ROWS:
    g = geocode(q)
    if g is None:  # fall back to tehsil/district level
        q2 = ", ".join(q.split(", ")[1:])
        g = geocode(q2)
        prec = "tehsil/district (village not found)"
    else:
        prec = "approx (district/town)" if "approx" in name else "village"
    out.append(dict(name=name, type=kind, query=q, capacity=cap, straw_t_yr=straw, precision=prec, lat=g[0], lon=g[1], matched=g[2]))
for name, lat, lon in COAL:
    out.append(dict(name=name, type="coal_cofiring", query="WRI GPPD", capacity=None, straw_t_yr=None, precision="plant", lat=lat, lon=lon, matched="WRI GPPD"))
df = pd.DataFrame(out)
df.to_csv(HERE / "buyers.csv", index=False)
gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.lon, df.lat), crs=4326).to_file(HERE / "buyers.geojson", driver="GeoJSON")
print(df[["name", "precision", "lat", "lon", "matched"]].to_string())
