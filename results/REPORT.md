# AnalogFinder — one-page report

**Question.** Which places on Earth physically resemble the Jezero, Gale and Meridiani landing sites, by how much, and where
does each analogy fail?

**Approach.** For each Mars target we build a fingerprint with five dimensions (topography at a matched 3.7 km / 18.5 km scale,
thermophysics, TIR mineralogy, aridity, thermal extremes) as a median and IQR over a 1°×1° window. We build the same fingerprint
for every 0.05° land cell on Earth from MODIS (MYD11C3, MCD43C3, MOD10CM), ASTER GED (AG1KM) and Copernicus GLO-90. Similarity
per dimension is s = exp(−d²/2), where d is the IQR-scaled distance to the target median. The overall score is a weighted mean that
requires all five dimensions. Thermal inertia and composition are compared as **within-planet percentiles**. The thermal-extremes
dimension in the weighted score is **Earth-relative** (how Mars-like a cell is compared with the rest of Earth's land); the absolute
comparison, which no Earth cell can pass, is shown separately as the "Mars gap". Details and limitations: `web/methods.html`.

## Artifact found and fixed (round 2)

**Symptom.** In round 1, 17 of Gale's 25 top candidate regions lay above 70°N (Greenland, Arctic Canada); Jezero had 1, Meridiani 0.
The Gale #1 card showed an Earth diurnal range of 3 K together with an apparent-thermal-inertia percentile of 0.98, and a TIR spectral
contrast of 0.014.

**Cause.** Apparent thermal inertia is (1 − albedo) / (LST_day − LST_night). Where the day–night contrast collapses (polar day and
polar night, near-coincident Aqua overpasses) it diverges, so the "bedrock-like" high ATI at high latitudes was a division artifact, not rock.
Snow-contaminated ASTER pixels have nearly flat emissivity spectra, so the mineral indices there were noise. A second, structural problem
surfaced while fixing the first: when the mineral dimension was unavailable its 30% weight was silently handed to the other dimensions,
so cells with no mineral signal could outrank cells with a real (and imperfect) one.

**Fix** (all thresholds in `config.yaml`, tests in `tests/`):
1. ATI is computed per month and a month counts only if LST_day − LST_night ≥ 10 K and LST_day > 0 °C; at least 4 usable months are
   required, otherwise ATI = n/a (`annual_ati`).
2. Mineral composition is n/a where TIR spectral contrast < 0.012 or the mean ASTER observation count (AG1km `NumObs`) < 3. The 0.012
   threshold separates the snow-flattened population (median contrast 0.007 at 60–70°N) from mid-latitude bare ground (median 0.019–0.068);
   the full table is in `results/mask_stats.txt`.
3. Permanent snow/ice from MOD10CM v061 (every usable month ≥ 50% snow, ≥ 3 usable months) is masked in addition to the warmest-month-LST rule.
4. The overall score requires all five weighted dimensions; a cell missing one is n/a instead of a free pass. Candidates and the web
   app use the same rule.

**Effect.** Valid land fell from 94.6% to 80.4% of land area; the ATI 95th percentile fell from 0.146 to 0.075 (the artifact tail).
Latitude histogram of the 25 top candidate regions per target, before → after:

| Target | -90..-60 | -60..-30 | -30..0 | 0..30 | 30..60 | 60..70 | 70..90 | >70N |
|---|---|---|---|---|---|---|---|---|
| Jezero Crater (before) | 0 | 0 | 0 | 6 | 5 | 13 | 1 | **1** |
| Jezero Crater (after) | 0 | 0 | 2 | 22 | 1 | 0 | 0 | **0** |
| Gale Crater (before) | 0 | 6 | 0 | 0 | 1 | 1 | 17 | **17** |
| Gale Crater (after) | 0 | 0 | 0 | 7 | 18 | 0 | 0 | **0** |
| Meridiani Planum (before) | 0 | 2 | 9 | 13 | 1 | 0 | 0 | **0** |
| Meridiani Planum (after) | 0 | 1 | 3 | 20 | 1 | 0 | 0 | **0** |

Real polar deserts do not survive, and we say so: Haughton Crater (75.4°N) and Svalbard are now n/a because fewer than four months have a
≥ 10 K day–night LST contrast there; Iceland's site cell is n/a for the same reason while cells within ±0.125° still score at the 92–97th
percentile. Haughton-like analogs would need a thermophysical proxy that does not depend on the diurnal cycle.


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


## Robustness to the weights

- Jezero Crater: top-25 regions stay candidates in a median 96% of weight samples (25 of 25 at ≥80%); they stay inside the top 25 in a median 44%.
- Gale Crater: top-25 regions stay candidates in a median 98% of weight samples (25 of 25 at ≥80%); they stay inside the top 25 in a median 42%.
- Meridiani Planum: top-25 regions stay candidates in a median 90% of weight samples (24 of 25 at ≥80%); they stay inside the top 25 in a median 38%.

Regions that only emerge under other weights are not re-labelled (the pool is the ~3,000 regions above threshold at the default weights).

## Blind predictions, checked with an instrument the score never used

Hematite, carbonates, clays and sulfates were not inputs to the global score (ASTER's five TIR bands cannot resolve them). After ranking,
each top-3 under-studied candidate was checked with EMIT L2B mineral identifications (60 m) and compared with 20,000 random pixels from
the same EMIT scenes. ✅ = present in ≥ 5% of candidate pixels and enriched ≥ 1.5× over the scene baseline. Everything is reported, including misses.

**Jezero Crater (Perseverance)** — expected from the Mars target: Calcite, Dolomite, Montmorillonite

| Candidate | Calcite (box vs scene baseline) | Dolomite (box vs scene baseline) | Montmorillonite (box vs scene baseline) |
|---|---|---|---|
| #1 Peru (-6.0, -81.1) | ✅ 25.4% vs 13.3% | ❌ 0.8% vs 0.6% | ❌ 0.1% vs 0.1% |
| #2 Saudi Arabia (20.5, 40.6) | ❌ 29.9% vs 37.8% | ❌ 1.2% vs 1.1% | ❌ 0.8% vs 1.6% |
| #3 Oman (24.5, 56.5) | ❌ 10.6% vs 33.5% | ❌ 9.6% vs 20.5% | ❌ 0.1% vs 0.3% |

**Gale Crater (Curiosity)** — expected from the Mars target: Montmorillonite, Kaolinite, Gypsum

| Candidate | Montmorillonite (box vs scene baseline) | Kaolinite (box vs scene baseline) | Gypsum (box vs scene baseline) |
|---|---|---|---|
| #1 China (37.6, 75.1) | ❌ 0.4% vs 0.6% | ❌ 1.5% vs 0.4% | ❌ 0.0% vs 0.0% |
| #2 China (40.5, 76.7) | ✅ 20.4% vs 4.7% | ❌ 0.0% vs 0.0% | ❌ 0.0% vs 0.0% |
| #3 Nepal (27.8, 86.9) | ❌ 0.7% vs 2.3% | ❌ 0.1% vs 2.3% | ❌ 0.0% vs 0.0% |

**Meridiani Planum (Opportunity)** — expected from the Mars target: Hematite, Gypsum

| Candidate | Hematite (box vs scene baseline) | Gypsum (box vs scene baseline) |
|---|---|---|
| #1 Australia (-24.2, 113.5) | ✅ 46.6% vs 19.8% | ✅ 13.1% vs 3.1% |
| #2 Kenya (2.6, 37.8) | ✅ 88.5% vs 41.0% | ❌ 0.0% vs 0.0% |
| #3 Oman (23.8, 57.3) | ❌ 46.6% vs 33.1% | ❌ 0.2% vs 0.2% |


## Case studies (Analog Cards)

- [Analog Card — China (37.63, 75.08) vs Gale Crater (Curiosity)](analog_cards/gale_1.md)
- [Analog Card — China (40.53, 76.66) vs Gale Crater (Curiosity)](analog_cards/gale_2.md)
- [Analog Card — Nepal (27.84, 86.90) vs Gale Crater (Curiosity)](analog_cards/gale_3.md)
- [Analog Card — Peru (-5.97, -81.09) vs Jezero Crater (Perseverance)](analog_cards/jezero_1.md)
- [Analog Card — Saudi Arabia (20.50, 40.58) vs Jezero Crater (Perseverance)](analog_cards/jezero_2.md)
- [Analog Card — Oman (24.55, 56.51) vs Jezero Crater (Perseverance)](analog_cards/jezero_3.md)
- [Analog Card — Australia (-24.19, 113.47) vs Meridiani Planum (Opportunity)](analog_cards/meridiani_1.md)
- [Analog Card — Kenya (2.64, 37.79) vs Meridiani Planum (Opportunity)](analog_cards/meridiani_2.md)
- [Analog Card — Oman (23.78, 57.26) vs Meridiani Planum (Opportunity)](analog_cards/meridiani_3.md)

## What works, what is approximate

- **Works:** both TIR-based dimensions use the same kind of measurement on both planets; the artifact above was caught by a latitude
  self-check and the fix is tested; the blind set and null test are reported as-is; the 1-D Mars thermal model reproduces Opportunity's
  measured 295 K peak (295.3 K).
- **Approximate:** the percentile bridges for thermal inertia and composition; Earth LST sampled only at the Aqua overpasses; Mars thermal
  extremes from a simple model; a 3.7 km topographic baseline set by the supplied MOLA grid; ASTER GED swath striping; MCD43C3 sampled on 24 days.
- **Meridiani remains the hardest target** for the classic analog list, which was assembled for other sites; its blind-set sites (Pilbara,
  Lonar, Dallol, Kilbourne Hole) score higher than its development set.
- **Choices made after seeing validation (round 1, disclosed):** bare-ground reference for Earth ATI percentiles; minimum scales 0.05 → 0.10
  (percentile features) and 0.02 → 0.05 (albedo); Copernicus GLO-90 instead of the gap-ridden ASTER GDEM. Round 2 changed masks and the
  missing-dimension rule, not tolerances; validation numbers were not tuned on.
- **Always breaks:** gravity (0.38 g vs 1 g), pressure (~6 vs ~1013 mbar), radiation, temperature (the Mars gap).
