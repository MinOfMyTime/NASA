# AnalogFinder — one-page report

**Question.** Which places on Earth physically resemble the Jezero, Gale and Meridiani landing sites, by how much, and where
does each analogy fail?

**Approach.** For each Mars target we build a fingerprint with five dimensions (topography at a matched 3.7 km / 18.5 km scale,
thermophysics, TIR mineralogy, aridity, thermal extremes) as a median and IQR over a 1°×1° window. We build the same fingerprint
for every 0.05° land cell on Earth from MODIS (MYD11C3, MCD43C3), ASTER GED (AG1KM) and Copernicus GLO-90. Similarity
per dimension is s = exp(−d²/2), where d is the IQR-scaled distance to the target median. The overall score is a weighted mean.
Thermal inertia and composition are compared as **within-planet percentiles**, because Earth's apparent thermal inertia and ASTER
band ratios are not in Mars units. Details and limitations: `web/methods.html`.

## Results

Validation: percentile of each known analog among all valid 0.05° land cells (cell at the site coordinate; best cell within ±0.125°).

| Known analog | Jezero Crater | Gale Crater | Meridiani Planum |
|---|---|---|---|
| Atacama Desert | 76 (best nearby 92) | 96 (best nearby 97) | 30 (best nearby 55) |
| Mojave / Death Valley | 48 (best nearby 91) | 98 (best nearby 99) | 47 (best nearby 79) |
| Utah MDRS | 33 (best nearby 79) | 68 (best nearby 88) | 36 (best nearby 63) |
| Iceland highlands | 86 (best nearby 98) | 82 (best nearby 92) | 55 (best nearby 79) |
| Mauna Kea flanks | 8 (best nearby 82) | 56 (best nearby 93) | 29 (best nearby 91) |
| Rio Tinto | 92 (best nearby 93) | 83 (best nearby 84) | 88 (best nearby 88) |
| Haughton Crater | 68 (best nearby 97) | 92 (best nearby 99) | 38 (best nearby 59) |

### Jezero Crater (Perseverance) — top 5 candidate regions (all within the top 0.5% of land)

Target fingerprint (1°×1° window, medians): TES TI 289 J m⁻² K⁻¹ s⁻½ (Mars percentile 0.82), albedo 0.175, relief 478 m, slope 1.5°, composition mafic/silica/fines 0.49/0.27/0.11, modelled surface T 181–288 K.

| # | Location (lat, lon) | Country | Area km² | Overall | Topo | Thermo | Mineral | Aridity | Extremes | Near known analog |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 19.17, 41.52 | Saudi Arabia | 17,253 | **0.81** | 0.88 | 0.90 | 0.81 | 0.80 | 0.00 | — |
| 2 | 65.04, 162.34 | Russia | 1,004 | **0.82** | 0.91 | 0.84 | 0.81 | 0.90 | 0.00 | — |
| 3 | 54.49, 158.42 | Russia | 126 | **0.83** | 0.94 | 0.93 | 0.85 | 0.74 | 0.00 | — |
| 4 | 62.89, 157.79 | Russia | 99 | **0.83** | 0.92 | 0.86 | 0.81 | 0.96 | 0.00 | — |
| 5 | 10.40, 11.47 | Nigeria | 122 | **0.83** | 0.94 | 0.93 | 0.86 | 0.65 | 0.00 | — |

### Gale Crater (Curiosity) — top 5 candidate regions (all within the top 0.5% of land)

Target fingerprint (1°×1° window, medians): TES TI 317 J m⁻² K⁻¹ s⁻½ (Mars percentile 0.89), albedo 0.223, relief 1531 m, slope 5.2°, composition mafic/silica/fines 0.38/0.39/0.18, modelled surface T 178–292 K.

| # | Location (lat, lon) | Country | Area km² | Overall | Topo | Thermo | Mineral | Aridity | Extremes | Near known analog |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 76.41, -25.19 | Greenland | 189 | **0.88** | 0.92 | 0.88 | 0.94 | 1.00 | 0.00 | — |
| 2 | 73.22, -81.30 | Canada | 36 | **0.90** | 0.93 | 0.95 | 0.95 | 0.95 | 0.00 | — |
| 3 | 82.29, -41.85 | Greenland | 41 | **0.89** | 0.95 | 0.95 | 0.90 | 1.00 | 0.00 | — |
| 4 | 76.72, -25.20 | Greenland | 71 | **0.89** | 0.94 | 0.95 | 0.89 | 1.00 | 0.00 | — |
| 5 | -51.34, -73.38 | Chile | 77 | **0.88** | 0.99 | 0.95 | 0.95 | 0.75 | 0.00 | — |

### Meridiani Planum (Opportunity) — top 5 candidate regions (all within the top 0.5% of land)

Target fingerprint (1°×1° window, medians): TES TI 223 J m⁻² K⁻¹ s⁻½ (Mars percentile 0.52), albedo 0.133, relief 77 m, slope 0.3°, composition mafic/silica/fines 0.55/0.33/0.11, modelled surface T 175–306 K.

| # | Location (lat, lon) | Country | Area km² | Overall | Topo | Thermo | Mineral | Aridity | Extremes | Near known analog |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | -24.19, 113.48 | Australia | 620 | **0.73** | 0.95 | 0.79 | 0.72 | 0.51 | 0.00 | — |
| 2 | -18.26, 123.94 | Australia | 323 | **0.70** | 0.95 | 0.91 | 0.68 | 0.20 | 0.00 | — |
| 3 | -24.15, 32.10 | Mozambique | 254 | **0.70** | 0.93 | 0.73 | 0.74 | 0.41 | 0.00 | — |
| 4 | -31.50, 130.00 | Australia | 343 | **0.69** | 0.83 | 0.81 | 0.80 | 0.29 | 0.00 | — |
| 5 | 44.78, 53.48 | Kazakhstan | 88 | **0.70** | 0.68 | 0.62 | 0.76 | 0.99 | 0.00 | — |


## Case studies (Analog Cards)

Top-3 under-studied candidates per target, plus Gale #4–#5 (#5, Patagonia, is the best Gale candidate inside EMIT coverage): ASTER GED 100 m, EMIT L2B minerals
(60 m) where coverage exists (EMIT is limited to roughly 54°S–54°N), and Copernicus GLO-30 slopes.

- [Analog Card — Greenland (76.41, -25.19) vs Gale Crater (Curiosity)](analog_cards/gale_1.md)
- [Analog Card — Canada (73.22, -81.30) vs Gale Crater (Curiosity)](analog_cards/gale_2.md)
- [Analog Card — Greenland (82.29, -41.85) vs Gale Crater (Curiosity)](analog_cards/gale_3.md)
- [Analog Card — Greenland (76.72, -25.20) vs Gale Crater (Curiosity)](analog_cards/gale_4.md)
- [Analog Card — Chile (-51.34, -73.38) vs Gale Crater (Curiosity)](analog_cards/gale_5.md)
- [Analog Card — Saudi Arabia (19.17, 41.52) vs Jezero Crater (Perseverance)](analog_cards/jezero_1.md)
- [Analog Card — Russia (65.04, 162.34) vs Jezero Crater (Perseverance)](analog_cards/jezero_2.md)
- [Analog Card — Russia (54.49, 158.42) vs Jezero Crater (Perseverance)](analog_cards/jezero_3.md)
- [Analog Card — Australia (-24.19, 113.48) vs Meridiani Planum (Opportunity)](analog_cards/meridiani_1.md)
- [Analog Card — Australia (-18.26, 123.94) vs Meridiani Planum (Opportunity)](analog_cards/meridiani_2.md)
- [Analog Card — Mozambique (-24.15, 32.10) vs Meridiani Planum (Opportunity)](analog_cards/meridiani_3.md)

**Independent check.** Hematite was *not* an input to the global score, because ASTER TIR cannot see it. Even so, EMIT finds hematite in 47% of
pixels at the #1 Meridiani candidate (Western Australia, 24.2°S 113.5°E) and in 70% at #2 (NW Australia, 18.3°S 123.9°E). Hematite is Meridiani's
signature mineral, so this is out-of-sample support for the method. At the #1 Jezero candidate (SW Saudi Arabia) EMIT's most common
mineral is calcite (46% of pixels); Jezero's delta is notable for carbonates.

## What works, what is approximate

- **Works:** both TIR-based dimensions use the same kind of measurement on both planets. Gale's classic analogs rank high (Mojave 99th,
  Atacama 96th, Haughton 92nd percentile), and for Jezero every listed analog has a cell at or above the 79th percentile within ±0.125°.
  The 1-D Mars thermal model reproduces Opportunity's measured 295 K peak (295.3 K).
- **Approximate:** the percentile bridges for thermal inertia and composition; Earth LST sampled only at the Aqua overpasses; Mars
  thermal extremes from a simple model; a 3.7 km topographic baseline set by the supplied MOLA grid; ASTER GED swath striping at high
  northern latitudes; MCD43C3 sampled on 24 days (each a 16-day retrieval window).
- **Meridiani is the hardest target.** It is a flat, dust-free, mafic, hematite-bearing plain. The listed analogs (chosen for other
  sites) mostly rank in the 29th–55th percentile; Rio Tinto, its classic chemical analog, ranks 88th. Its best candidates are flat,
  iron-oxide-rich plains in Australia, plus Mozambique and Kazakhstan.
- **Choices made after seeing the validation (disclosed):** Earth thermal-inertia percentiles use bare ground (NDVI < 0.3) as the
  reference population; the minimum scale for percentile features was raised from 0.05 to 0.10 and for albedo from 0.02 to 0.05;
  ASTER GED elevation was replaced by Copernicus GLO-90 because the GDEM layer is all-fill in 4,944 of 24,873 tiles (including Hawaii).
- **Always breaks:** gravity (0.38 g vs 1 g), pressure (~6 vs ~1013 mbar), radiation, temperature. The thermal-extremes score is
  0.00 for every Earth cell we checked, which is shown as a "breaks down" dimension.
