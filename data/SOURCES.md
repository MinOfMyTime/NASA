# Data sources

All raw files live in `data/raw/` (git-ignored). Access dates are UTC.

## Mars

| File | Product | Source / DOI | Version | Access |
|---|---|---|---|---|
| `megt90n000eb.img` | MOLA MEGDR median topography, 16 ppd (1/16 deg, 3.705 km), MSB int16 m, line 1 = 90N, sample 1 = 0E, east-positive | PDS Geosciences Node, `MGS-M-MOLA-5-MEGDR-L3-V1.0` — https://pds-geosciences.wustl.edu/mgs/mgs-m-mola-5-megdr-l3-v1/mgsl_300x/meg016/ (Smith et al. 2003, doi:10.17189/1519460) | PRODUCT_VERSION_ID 2.0 | supplied by user 2026-11-14 |
| `Mars_MGS_TES_Albedo_mosaic_global_7410m.tif` | TES bolometric albedo, 8 ppd, simple cylindrical, -180..180E | USGS Astrogeology — https://planetarymaps.usgs.gov/mosaic/Mars_MGS_TES_Albedo_mosaic_global_7410m.tif (Christensen et al. 2001, doi:10.1029/2000JE001370) | — | supplied by user 2026-11-14 |
| `nmap2003.bin`, `nmsk2003.bin` | TES nightside thermal inertia, 20 ppd, big-endian int16 (24–800 J m⁻² K⁻¹ s⁻½, 0 = no data) + interpolation mask. Header: 7200-byte text, latitudes -90→90 (south first), **longitude 0–360 WEST** | Putzig & Mellon (2005) Icarus 173, 325, doi:10.1016/j.icarus.2004.08.017 — http://www.mars.asu.edu/data/ / se.psi.edu/inertia/ti2003.html | 2003-06-13 map | supplied by user 2026-11-14 |
| `tes_minerals/tes_*.vic` | TES mineral maps (plagioclase, high/low-Ca pyroxene, olivine, quartz, sheet silicates/high-Si glass, surface dust, hematite), 4 ppd, VICAR REAL little-endian, -1 = no data | ASU Mars Space Flight Facility — https://www.mars.asu.edu/data/ (Bandfield 2002, JGR 107, doi:10.1029/2001JE001510) | VICAR conversion 2026-11-14 (label `IOMEDLEY`) | supplied by user 2026-11-14 |

Orientation checks (`scripts/inspect_mars.py`): TI header orientation gives corr(TI, albedo) = −0.42 (alternatives ≤ −0.12);
the VICAR labels carry no map geometry, so north-up / 0°E-origin was established by corr(dust, albedo) = +0.62 (best of 8
hypotheses) and confirmed by TES hematite at Meridiani (window median 0.125 > global 99th percentile 0.068).

## Earth (NASA Earthdata / LP DAAC, accessed 2026-11-14 with an Earthdata Login token)

| Dir | Product | DOI | Version | Use |
|---|---|---|---|---|
| `MYD11C3/` | Aqua MODIS LST/Emissivity Monthly L3 Global 0.05 deg CMG, 12 granules for 2023 | doi:10.5067/MODIS/MYD11C3.061 | 061 | LST_Day_CMG, LST_Night_CMG (uint16 × 0.02 K, fill 0), QC_Day/QC_Night |
| `MCD43C3/` | MODIS BRDF/Albedo Daily L3 0.05 deg CMG, days 1 and 16 of every month of 2023 (24 granules; each is a 16-day moving-window retrieval, so the samples cover the whole year) | doi:10.5067/MODIS/MCD43C3.061 | 061 | Albedo_WSA_shortwave (int16 × 0.001, fill 32767), Albedo_Quality ≤ 3, Percent_Snow = 0 |
| `AG1KM/` | ASTER Global Emissivity Dataset 1 km, all 1×1 deg tiles (24,873) | doi:10.5067/ASTER/AG1KM.003 | 003 | Emissivity/Mean bands 10–14 (× 0.001), NDVI/Mean (× 0.01), ASTER GDEM (m), Land Water Map (0 land, 1 water) |
| `AG100/` | ASTER GED 100 m tiles (case studies only) | doi:10.5067/ASTER/AG100.003 | 003 | zoom-ins |
| `EMITL2BMIN/` | EMIT L2B Estimated Mineral Identification and Band Depth, 60 m (case studies only) | doi:10.5067/EMIT/EMITL2BMIN.001 | 001 | zoom-ins |
| (no local file) `data/processed/earth_dem_30ppd.npy` | Copernicus DEM GLO-90 (public AWS bucket `copernicus-dem-90m`, 26,475 tiles, 1/4 COG overview area-averaged to 1/30 deg by `src/download_dem.py`) — https://copernicus-dem-90m.s3.amazonaws.com/ , doi:10.5270/ESA-c5d3d65 | 2021_1 | Earth topography (replaces ASTER GED GDEM, which is all-fill in 4,944 tiles) |
| (streamed) | Copernicus DEM GLO-30 (public AWS bucket `copernicus-dem-30m`) | doi:10.5270/ESA-c5d3d65 | 2021_1 | case-study slopes |
| `emit_mineral_grouping_matrix_20230503.csv` | EMIT L2B mineral grouping matrix (maps mineral IDs to the 10 core minerals) — https://github.com/emit-sds/emit-sds-l2b/blob/develop/data/mineral_grouping_matrix_20230503.csv | — | 20230503 | case studies |
| `ne_50m_admin_0_countries.zip` | Natural Earth admin-0 countries 1:50m (public domain) — https://naciscdn.org/naturalearth/50m/cultural/ne_50m_admin_0_countries.zip | — | 5.x | country labels for candidates |

## Literature used for checks (not as data inputs)

- Spanovich et al. (2006) Icarus 180, 314, doi:10.1016/j.icarus.2005.09.014 — "Opportunity experienced a maximum daily surface
  temperature of 295 K … at the beginning of the mission (Ls = 339°)". Our 1-D model gives 295.3 K for the Meridiani window
  (`scripts/check_thermal_model.py`).
- Ninomiya, Fu & Cudahy (2005) Remote Sens. Environ. 99, 127, doi:10.1016/j.rse.2005.06.009 — ASTER TIR quartz and mafic indices.
