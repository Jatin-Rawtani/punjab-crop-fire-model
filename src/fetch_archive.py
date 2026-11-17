"""Download NASA FIRMS yearly country archives for India (public, no key) and keep only
detections inside a Punjab + Haryana box. Writes raw_archive/<sensor>_<year>_PBHR.csv.

Sensors: VIIRS S-NPP 375 m (2012-2025, one consistent sensor for trends)
         MODIS 1 km     (2001-2025, the long record)
Source:  https://firms.modaps.eosdis.nasa.gov/data/country/  (NASA FIRMS, open data)
Run:     python src/fetch_archive.py
"""
import io, pathlib, sys, urllib.request
import pandas as pd

OUT = pathlib.Path(__file__).resolve().parent.parent / "data" / "raw_archive"
OUT.mkdir(parents=True, exist_ok=True)
BOX = dict(lon=(73.8, 77.7), lat=(27.6, 32.6))  # Punjab + Haryana
JOBS = [("viirs-snpp", y) for y in range(2012, 2026)] + [("modis", y) for y in range(2001, 2026)]

for sensor, year in JOBS:
    dest = OUT / f"{sensor}_{year}_PBHR.csv"
    if dest.exists():
        continue
    url = f"https://firms.modaps.eosdis.nasa.gov/data/country/{sensor}/{year}/{sensor}_{year}_India.csv"
    try:
        raw = urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "punjab-crop-fire-model"}), timeout=300).read()
    except Exception as e:
        print("FAIL", sensor, year, e, flush=True)
        continue
    df = pd.read_csv(io.BytesIO(raw))
    m = df.longitude.between(*BOX["lon"]) & df.latitude.between(*BOX["lat"])
    df[m].to_csv(dest, index=False)
    print(sensor, year, "india", len(df), "box", int(m.sum()), flush=True)
print("done")
