# AnalogFinder — a cross-planet fingerprint search engine for Mars analogs on Earth

NASA Space Apps Challenge 2026 · *Find and characterize new terrestrial analogs for Moon/Mars landing sites and bases.*

**Live app: https://minofmytime.github.io/NASA/** (static; deployed from `web/` by GitHub Actions on every push to `main`)

Most analog work compares photos. **AnalogFinder compares physics.** For a Mars target (Jezero, Gale, Meridiani) we
build a five-dimension fingerprint from orbital Mars data. We then build the *same* fingerprint for every 0.05° land cell on Earth
from NASA Earth data and rank the planet, per dimension and overall. Both planets have been mapped in the
**thermal infrared** (TES on Mars, MODIS and ASTER on Earth), so the thermophysics and mineralogy comparisons rest on the
same kind of measurement, not on appearance.

Every result states **what matches**, **how well** (0–1 per dimension) and **where the analogy breaks** (gravity, pressure,
radiation, temperature, biology).

| Dimension | Mars | Earth |
|---|---|---|
| Topography (3.7 km baseline, 18.5 km windows) | MOLA MEGDR 1/16° | Copernicus GLO-90 → 1/30° |
| Thermophysics | TES thermal inertia + albedo | MODIS apparent thermal inertia (MYD11C3 day/night LST + MCD43C3 albedo) + albedo |
| Mineralogy (TIR) | TES mineral maps → mafic / silica / dust | ASTER GED bands 10–14 → mafic index, quartz index, spectral contrast |
| Aridity | NDVI ≡ 0 | ASTER GED NDVI |
| Thermal extremes | 1-D thermal model from TES TI + albedo | MODIS monthly LST |

See [`web/methods.html`](web/methods.html) for formulas and limitations, and [`data/SOURCES.md`](data/SOURCES.md) for every file's URL, DOI, version and access date.

## Results

Validation: percentile of each site among all valid 0.05° land cells (cell at the site coordinate; best cell within ±0.125°).

**Development set** (used while building the method):

| Site | Jezero Crater | Gale Crater | Meridiani Planum |
|---|---|---|---|
| Atacama Desert | 88 (best nearby 98) | 98 (best nearby 99) | 36 (best nearby 69) |
| Mojave / Death Valley | 75 (best nearby 96) | 100 (best nearby 100) | 54 (best nearby 81) |
| Utah MDRS | 48 (best nearby 91) | 80 (best nearby 94) | 49 (best nearby 73) |
| Iceland highlands | nan (best nearby 96) | nan (best nearby 92) | nan (best nearby 62) |
| Mauna Kea flanks | 19 (best nearby 91) | 69 (best nearby 95) | 42 (best nearby 93) |
| Rio Tinto | 97 (best nearby 97) | 89 (best nearby 89) | 70 (best nearby 71) |
| Haughton Crater | nan (best nearby nan) | nan (best nearby nan) | nan (best nearby nan) |

**Blind set** (published analogs scored once with the frozen settings; references in `results/validation.md`):

| Site | Jezero Crater | Gale Crater | Meridiani Planum |
|---|---|---|---|
| Pilbara, North Pole Dome (Dresser Fm.) | 88 (best nearby 93) | 90 (best nearby 93) | 98 (best nearby 99) |
| McMurdo Dry Valleys (Taylor Valley) | nan (best nearby nan) | nan (best nearby nan) | nan (best nearby nan) |
| Lonar Crater | nan (best nearby nan) | nan (best nearby nan) | nan (best nearby nan) |
| Dallol / Danakil | nan (best nearby 100) | nan (best nearby 100) | nan (best nearby 99) |
| Svalbard, Bockfjord volcanic complex (Sverrefjell) | nan (best nearby nan) | nan (best nearby nan) | nan (best nearby nan) |
| Kilbourne Hole maar | 89 (best nearby 95) | 84 (best nearby 93) | 91 (best nearby 97) |

**Null test** (1,000 random valid land cells, seed 0):

- Jezero Crater: median percentile 48; 24% of random cells reach the 75th percentile, 11% the 90th; histogram [0-25, 25-50, 50-75, 75-90, 90-100] = [262, 266, 231, 133, 108]
- Gale Crater: median percentile 51; 24% of random cells reach the 75th percentile, 9% the 90th; histogram [0-25, 25-50, 50-75, 75-90, 90-100] = [255, 229, 275, 149, 92]
- Meridiani Planum: median percentile 49; 25% of random cells reach the 75th percentile, 11% the 90th; histogram [0-25, 25-50, 50-75, 75-90, 90-100] = [262, 255, 236, 139, 108]

### Jezero Crater (Perseverance) — top 5 candidate regions (top 0.5% of land, ranked by mean overall score)

Target fingerprint (1°×1° window, medians): TES TI 289 J m⁻² K⁻¹ s⁻½ (Mars percentile 0.82), albedo 0.175, relief 478 m, slope 1.5°, composition mafic/silica/fines 0.49/0.27/0.11, modelled surface T 181–288 K.

| # | Location (lat, lon) | Country | Area km² | Overall | Topo | Thermo | Mineral | Aridity | Extremes (rel.) | Robust: candidate / top-25 |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | -5.97, -81.09 | Peru | 154 | **0.84** | 0.92 | 0.92 | 0.72 | 0.98 | 0.32 | 99% / 73% |
| 2 | 20.50, 40.58 | Saudi Arabia | 5,646 | **0.84** | 0.87 | 0.89 | 0.84 | 0.86 | 0.28 | 99% / 73% |
| 3 | 24.55, 56.51 | Oman | 112 | **0.82** | 0.70 | 0.94 | 0.83 | 0.97 | 0.33 | 92% / 53% |
| 4 | 14.47, 40.75 | Eritrea | 120 | **0.82** | 0.90 | 0.89 | 0.72 | 0.93 | 0.32 | 98% / 53% |
| 5 | -13.32, 12.70 | Angola | 120 | **0.82** | 0.95 | 0.94 | 0.73 | 0.72 | 0.38 | 98% / 47% |

Robustness = share of 500 Dirichlet-sampled weight vectors (centred on the defaults) in which the region still scores above the top-0.5% threshold / still ranks in the top 25. Extremes (rel.) is Earth-relative; the absolute "Mars gap" is reported separately.

### Gale Crater (Curiosity) — top 5 candidate regions (top 0.5% of land, ranked by mean overall score)

Target fingerprint (1°×1° window, medians): TES TI 317 J m⁻² K⁻¹ s⁻½ (Mars percentile 0.89), albedo 0.223, relief 1531 m, slope 5.2°, composition mafic/silica/fines 0.38/0.39/0.18, modelled surface T 178–292 K.

| # | Location (lat, lon) | Country | Area km² | Overall | Topo | Thermo | Mineral | Aridity | Extremes (rel.) | Robust: candidate / top-25 |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 58.33, -155.15 | United States of America | 65 | **0.91** | 0.98 | 0.95 | 0.94 | 0.83 | 0.40 | 100% / 98% |
| 2 | 37.63, 75.08 | China | 171 | **0.88** | 0.96 | 0.71 | 0.97 | 0.86 | 0.72 | 100% / 88% |
| 3 | 40.53, 76.66 | China | 141 | **0.86** | 0.88 | 0.84 | 0.95 | 0.74 | 0.72 | 100% / 70% |
| 4 | 27.84, 86.90 | Nepal | 191 | **0.86** | 0.92 | 0.79 | 0.97 | 0.76 | 0.57 | 100% / 68% |
| 5 | 43.04, 94.34 | China | 226 | **0.86** | 0.89 | 0.89 | 0.97 | 0.63 | 0.65 | 98% / 68% |

Robustness = share of 500 Dirichlet-sampled weight vectors (centred on the defaults) in which the region still scores above the top-0.5% threshold / still ranks in the top 25. Extremes (rel.) is Earth-relative; the absolute "Mars gap" is reported separately.

### Meridiani Planum (Opportunity) — top 5 candidate regions (top 0.5% of land, ranked by mean overall score)

Target fingerprint (1°×1° window, medians): TES TI 223 J m⁻² K⁻¹ s⁻½ (Mars percentile 0.52), albedo 0.133, relief 77 m, slope 0.3°, composition mafic/silica/fines 0.55/0.33/0.11, modelled surface T 175–306 K.

| # | Location (lat, lon) | Country | Area km² | Overall | Topo | Thermo | Mineral | Aridity | Extremes (rel.) | Robust: candidate / top-25 |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | -24.19, 113.47 | Australia | 508 | **0.75** | 0.95 | 0.63 | 0.87 | 0.49 | 0.46 | 96% / 80% |
| 2 | 2.64, 37.79 | Kenya | 185 | **0.74** | 0.87 | 0.61 | 0.71 | 0.96 | 0.30 | 95% / 83% |
| 3 | 23.78, 57.26 | Oman | 368 | **0.72** | 0.80 | 0.82 | 0.49 | 0.98 | 0.45 | 98% / 66% |
| 4 | -31.49, 130.01 | Australia | 395 | **0.71** | 0.82 | 0.87 | 0.72 | 0.30 | 0.62 | 90% / 64% |
| 5 | 44.78, 53.48 | Kazakhstan | 88 | **0.71** | 0.68 | 0.59 | 0.71 | 0.99 | 0.61 | 88% / 57% |

Robustness = share of 500 Dirichlet-sampled weight vectors (centred on the defaults) in which the region still scores above the top-0.5% threshold / still ranks in the top 25. Extremes (rel.) is Earth-relative; the absolute "Mars gap" is reported separately.

Full report (incl. the artifact found and fixed, robustness and blind predictions): [`results/REPORT.md`](results/REPORT.md) · validation: [`results/validation.md`](results/validation.md) · candidates: `results/candidates_<target>.csv` · cards: [`results/analog_cards/`](results/analog_cards/)

![AnalogFinder web app](results/screenshot_app.png)

## Run from scratch

```bash
python3 -m venv .venv && . .venv/bin/activate && pip install -r requirements.txt
# Earthdata access: put an EDL token in ~/.earthdata_token (or export EARTHDATA_TOKEN),
# or EARTHDATA_USERNAME/EARTHDATA_PASSWORD, or ~/.netrc.
# Mars inputs go in data/raw/mars/ (see data/SOURCES.md for the URLs):
#   megt90n000eb.img, Mars_MGS_TES_Albedo_mosaic_global_7410m.tif, nmap2003.bin, nmsk2003.bin, tes_minerals/*.vic

python scripts/inspect_mars.py            # MVP-0: Mars shapes/ranges + orientation checks + quicklooks
python src/download_earth.py modis        # MYD11C3 x12 + MCD43C3 x24 (~5.6 GB)
python src/download_earth.py snow         # MOD10CM x12 (permanent snow/ice mask)
python src/download_earth.py ag1km        # ASTER GED 1 km, 24,873 tiles (~11 GB)
python src/download_dem.py                # Copernicus GLO-90 -> 1/30 deg (public AWS, overview reads)
python scripts/quicklook_earth.py         # MVP-0: Earth quicklooks
python src/features_mars.py               # Mars features + target fingerprints (+ thermal model)
python scripts/check_thermal_model.py     # model vs Opportunity's observed 295 K
python src/features_earth.py              # Earth features on the 0.05 deg grid
python src/similarity.py                  # per-dimension + overall similarity rasters
python src/validate.py                    # results/validation.md
python src/candidates.py                  # results/candidates_<target>.csv
python scripts/lat_hist.py after          # latitude self-check of the candidate lists
python src/robustness.py                  # 500 weight samples -> robustness badges (~25 min per target)
python src/case_study.py jezero 3         # Analog Cards (AG100 + EMIT + GLO-30); also gale / meridiani
python src/export_web.py                  # web/data/*
python scripts/make_report.py             # results/REPORT.md + README results tables
python -m http.server -d web 8000         # open http://localhost:8000
pytest -q                                 # unit tests
```

Runtime on 4 vCPU / 16 GB: downloads ~30 min, Earth features ~10 min, everything else a few minutes.

## Repository

```
config.yaml            targets (Mars lon east-positive 0-360), windows, weights, scales, masks
src/                   geo, io_mars, topo (identical on both planets), stats, features_*, similarity, validate,
                       candidates, export_web, case_study, download_*, mars_thermal_model
scripts/               inspection, quicklooks and debugging aids
results/               validation.md, candidates_<target>.csv, analog_cards/, quicklooks/, REPORT.md
web/                   static single-page app (Leaflet), deployable to GitHub Pages as-is
tests/                 longitude conventions, orientation flips, scale factors/QC bits, ATI, percentiles, similarity
```

## Data credits

NASA LP DAAC (MODIS MYD11C3/MCD43C3 v061, ASTER GED v003, EMIT L2B v001); NASA PDS Geosciences Node (MOLA MEGDR);
USGS Astrogeology (TES albedo mosaic); Putzig & Mellon (TES thermal inertia); ASU Mars Space Flight Facility (TES mineral
maps); ESA Copernicus DEM GLO-90/GLO-30 (© DLR e.V. 2010-2014 and © Airbus Defence and Space GmbH 2014-2018, provided under
COPERNICUS by the European Union and ESA); Natural Earth; basemap © Esri, HERE, Garmin, © OpenStreetMap contributors.
