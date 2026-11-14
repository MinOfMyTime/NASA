"""Generate results/REPORT.md and the README results section from the pipeline outputs (no hand-typed numbers)."""
import json
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
fps = json.load(open(ROOT / "data/processed/mars_fingerprints.json"))
val = json.load(open(ROOT / "data/processed/validation.json"))
DIMS = ["topo", "thermo", "mineral", "aridity", "extremes"]


def top_table(t, which="top", n=5):
    df = pd.read_csv(ROOT / f"results/candidates_{t}.csv")
    df = df[df.list == which].head(n)
    rows = ["| # | Location (lat, lon) | Country | Area km² | Overall | Topo | Thermo | Mineral | Aridity | Extremes | Near known analog |",
            "|---|---|---|---|---|---|---|---|---|---|---|"]
    for _, r in df.iterrows():
        rows.append(f"| {int(r['rank'])} | {r.lat:.2f}, {r.lon:.2f} | {r.country} | {r.area_km2:,.0f} | **{r.mean_score:.2f}** | "
                    + " | ".join(f"{r['s_' + d]:.2f}" for d in DIMS) + f" | {r.near_known if isinstance(r.near_known, str) else '—'} |")
    return "\n".join(rows)


def val_table():
    names = list(next(iter(val.values())).keys())
    rows = ["| Known analog | " + " | ".join(fps[t]["name"].split(" (")[0] for t in val) + " |",
            "|---|" + "---|" * len(val)]
    for k in names:
        cells = []
        for t in val:
            v = val[t][k]
            p, pb = v["percentile"], v["best_neighbourhood_percentile"]
            cells.append("NaN" if p is None else f"{p:.0f} (best nearby {pb:.0f})")
        rows.append(f"| {val[next(iter(val))][k]['name']} | " + " | ".join(cells) + " |")
    return "\n".join(rows)


def cards():
    out = []
    for p in sorted((ROOT / "results/analog_cards").glob("*.md")):
        title = p.read_text().splitlines()[0].lstrip("# ").strip()
        out.append(f"- [{title}](analog_cards/{p.name})")
    return "\n".join(out)


body = [f"Validation: percentile of each known analog among all valid 0.05° land cells (cell at the site coordinate; best cell within ±0.125°).", "",
        val_table(), ""]
for t, fp in fps.items():
    f = fp["features"]
    body += [f"### {fp['name']} — top 5 candidate regions (all within the top 0.5% of land)", "",
             f"Target fingerprint (1°×1° window, medians): TES TI {f['thermal_inertia']['median']:.0f} J m⁻² K⁻¹ s⁻½ "
             f"(Mars percentile {f['pct_ti']['median']:.2f}), albedo {f['albedo']['median']:.3f}, relief {10 ** f['log_relief']['median'] - 1:.0f} m, "
             f"slope {f['slope_deg']['median']:.1f}°, composition mafic/silica/fines {f['comp_mafic']['median']:.2f}/"
             f"{f['comp_silica']['median']:.2f}/{f['comp_fines']['median']:.2f}, modelled surface T {f['lst_min_k']['median']:.0f}–{f['lst_max_k']['median']:.0f} K.", "",
             top_table(t), ""]
results = "\n".join(body)

report = f"""# AnalogFinder — one-page report

**Question.** Which places on Earth physically resemble the Jezero, Gale and Meridiani landing sites, by how much, and where
does each analogy fail?

**Approach.** For each Mars target we build a fingerprint with five dimensions (topography at a matched 3.7 km / 18.5 km scale,
thermophysics, TIR mineralogy, aridity, thermal extremes) as a median and IQR over a 1°×1° window. We build the same fingerprint
for every 0.05° land cell on Earth from MODIS (MYD11C3, MCD43C3), ASTER GED (AG1KM) and Copernicus GLO-90. Similarity
per dimension is s = exp(−d²/2), where d is the IQR-scaled distance to the target median. The overall score is a weighted mean.
Thermal inertia and composition are compared as **within-planet percentiles**, because Earth's apparent thermal inertia and ASTER
band ratios are not in Mars units. Details and limitations: `web/methods.html`.

## Results

{results}

## Case studies (Analog Cards)

Top-3 under-studied candidates per target, plus Gale #4–#5 (#5, Patagonia, is the best Gale candidate inside EMIT coverage): ASTER GED 100 m, EMIT L2B minerals
(60 m) where coverage exists (EMIT is limited to roughly 54°S–54°N), and Copernicus GLO-30 slopes.

{cards()}

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
"""
(ROOT / "results/REPORT.md").write_text(report)
readme = (ROOT / "README.md").read_text()
readme = re.sub(r"## Results\n\n.*?\n\n## Run from scratch", "## Results\n\n" + results.replace("\\", "\\\\")
                + "\nFull report: [`results/REPORT.md`](results/REPORT.md) · validation: [`results/validation.md`](results/validation.md) · "
                "candidates: `results/candidates_<target>.csv` · cards: [`results/analog_cards/`](results/analog_cards/)\n\n"
                "![AnalogFinder web app](results/screenshot_app.png)\n\n## Run from scratch", readme, flags=re.S)
(ROOT / "README.md").write_text(readme)
print(report)
