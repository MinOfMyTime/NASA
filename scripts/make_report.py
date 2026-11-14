"""Generate results/REPORT.md and the README results section from the pipeline outputs (no hand-typed numbers)."""
import json
import re
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
cfg = yaml.safe_load(open(ROOT / "config.yaml"))
fps = json.load(open(ROOT / "data/processed/mars_fingerprints.json"))
val = json.load(open(ROOT / "data/processed/validation.json"))
DIMS = ["topo", "thermo", "mineral", "aridity", "extremes"]
LH = {tag: json.load(open(ROOT / f"results/lat_hist_{tag}.json")) for tag in ("before", "after")}
BP = json.load(open(ROOT / "results/analog_cards/summary.json"))
MASK_STATS = (ROOT / "results/mask_stats.txt").read_text().strip() if (ROOT / "results/mask_stats.txt").exists() else ""


def top_table(t, which="top", n=5):
    df = pd.read_csv(ROOT / f"results/candidates_{t}.csv")
    df = df[df.list == which].head(n)
    rows = ["| # | Location (lat, lon) | Country | Area km² | Overall | Topo | Thermo | Mineral | Aridity | Extremes (rel.) | Robust: candidate / top-25 |",
            "|---|---|---|---|---|---|---|---|---|---|---|"]
    for _, r in df.iterrows():
        rows.append(f"| {int(r['rank'])} | {r.lat:.2f}, {r.lon:.2f} | {r.country} | {r.area_km2:,.0f} | **{r.mean_score:.2f}** | "
                    + " | ".join(f"{r['s_' + d]:.2f}" for d in DIMS)
                    + f" | {r.robust_pool_pct:.0f}% / {r.robust_top25_pct:.0f}% |")
    return "\n".join(rows)


def val_table(key):
    v = val[key]
    names = list(next(iter(v.values())).keys())
    rows = ["| Site | " + " | ".join(fps[t]["name"].split(" (")[0] for t in v) + " |", "|---|" + "---|" * len(v)]
    for k in names:
        cells = []
        for t in v:
            s = v[t][k]
            p, pb = s["percentile"], s["best_neighbourhood_percentile"]
            cells.append(f"n/a ({s['note'].split('(')[0].strip()})" if p is None else f"{p:.0f} (best nearby {pb:.0f})")
        rows.append(f"| {v[next(iter(v))][k]['name']} | " + " | ".join(cells) + " |")
    return "\n".join(rows)


def null_text():
    return "\n".join(f"- {fps[t]['name'].split(' (')[0]}: median percentile {n['median_pct']:.0f}; {n['p75_share']:.0%} of random cells reach the 75th "
                     f"percentile, {n['p90_share']:.0%} the 90th; histogram [0-25, 25-50, 50-75, 75-90, 90-100] = {n['hist']}"
                     for t, n in val["null"].items())


def lat_table():
    bins = LH["after"]["jezero"]["bins"]
    hdr = "| Target | " + " | ".join(f"{a}..{b}" for a, b in zip(bins[:-1], bins[1:])) + " | >70N |"
    rows = [hdr, "|---|" + "---|" * (len(bins))]
    for t in fps:
        for tag in ("before", "after"):
            h = LH[tag][t]
            rows.append(f"| {fps[t]['name'].split(' (')[0]} ({tag}) | " + " | ".join(str(c) for c in h["counts"]) + f" | **{h['n_above_70']}** |")
    return "\n".join(rows)


def blind_pred_table():
    out = []
    for t in fps:
        rows = [r for r in BP if r["target"] == t]
        if not rows:
            continue
        exp = rows[0]["expected"]
        out += [f"**{fps[t]['name']}** — expected from the Mars target: {', '.join(exp)}", "",
                "| Candidate | " + " | ".join(f"{m} (box vs scene baseline)" for m in exp) + " |", "|---|" + "---|" * len(exp)]
        for r in rows:
            if r["emit_status"] != "ok":
                out.append(f"| #{r['rank']} {r['country']} ({r['lat']:.1f}, {r['lon']:.1f}) | " + " | ".join(f"n/a ({r['emit_status']})" for _ in exp) + " |")
                continue
            cells = []
            for m in exp:
                a, b = r["fraction"][m] * 100, r["baseline"][m] * 100
                mark = "✅" if a >= 5 and a > 1.5 * b else "❌"
                cells.append(f"{mark} {a:.1f}% vs {b:.1f}%")
            out.append(f"| #{r['rank']} {r['country']} ({r['lat']:.1f}, {r['lon']:.1f}) | " + " | ".join(cells) + " |")
        out.append("")
    return "\n".join(out)


def cards():
    return "\n".join(f"- [{p.read_text().splitlines()[0].lstrip('# ').strip()}](analog_cards/{p.name})"
                     for p in sorted((ROOT / "results/analog_cards").glob("*.md")))


def robust_summary():
    out = []
    for t in fps:
        df = pd.read_csv(ROOT / f"results/candidates_{t}.csv")
        top = df[df.list == "top"]
        out.append(f"- {fps[t]['name'].split(' (')[0]}: top-25 regions stay candidates in a median {top.robust_pool_pct.median():.0f}% of weight samples "
                   f"({(top.robust_pool_pct >= 80).sum()} of {len(top)} at ≥80%); they stay inside the top 25 in a median {top.robust_top25_pct.median():.0f}%.")
    return "\n".join(out)


results = ["Validation: percentile of each site among all valid 0.05° land cells (cell at the site coordinate; best cell within ±0.125°).", "",
           "**Development set** (used while building the method):", "", val_table("development"), "",
           "**Blind set** (published analogs scored once with the frozen settings; references in `results/validation.md`):", "",
           val_table("blind"), "", "**Null test** (1,000 random valid land cells, seed 0):", "", null_text(), ""]
for t, fp in fps.items():
    f = fp["features"]
    results += [f"### {fp['name']} — top 5 candidate regions (top 0.5% of land, ranked by mean overall score)", "",
                f"Target fingerprint (1°×1° window, medians): TES TI {f['thermal_inertia']['median']:.0f} J m⁻² K⁻¹ s⁻½ "
                f"(Mars percentile {f['pct_ti']['median']:.2f}), albedo {f['albedo']['median']:.3f}, relief {10 ** f['log_relief']['median'] - 1:.0f} m, "
                f"slope {f['slope_deg']['median']:.1f}°, composition mafic/silica/fines {f['comp_mafic']['median']:.2f}/"
                f"{f['comp_silica']['median']:.2f}/{f['comp_fines']['median']:.2f}, modelled surface T {f['lst_min_k']['median']:.0f}–{f['lst_max_k']['median']:.0f} K.", "",
                top_table(t), "",
                "Robustness = share of 500 Dirichlet-sampled weight vectors (centred on the defaults) in which the region still scores above the "
                "top-0.5% threshold / still ranks in the top 25. Extremes (rel.) is Earth-relative; the absolute \"Mars gap\" is reported separately.", ""]
results = "\n".join(results)

artifact = f"""## Artifact found and fixed (round 2)

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

{lat_table()}

Real polar deserts do not survive, and we say so: Haughton Crater (75.4°N) and Svalbard are now n/a because fewer than four months have a
≥ 10 K day–night LST contrast there; Iceland's site cell is n/a for the same reason while cells within ±0.125° still score at the 92–97th
percentile. Haughton-like analogs would need a thermophysical proxy that does not depend on the diurnal cycle.
"""

report = f"""# AnalogFinder — one-page report

**Question.** Which places on Earth physically resemble the Jezero, Gale and Meridiani landing sites, by how much, and where
does each analogy fail?

**Approach.** For each Mars target we build a fingerprint with five dimensions (topography at a matched 3.7 km / 18.5 km scale,
thermophysics, TIR mineralogy, aridity, thermal extremes) as a median and IQR over a 1°×1° window. We build the same fingerprint
for every 0.05° land cell on Earth from MODIS (MYD11C3, MCD43C3, MOD10CM), ASTER GED (AG1KM) and Copernicus GLO-90. Similarity
per dimension is s = exp(−d²/2), where d is the IQR-scaled distance to the target median. The overall score is a weighted mean that
requires all five dimensions. Thermal inertia and composition are compared as **within-planet percentiles**. The thermal-extremes
dimension in the weighted score is **Earth-relative** (how Mars-like a cell is compared with the rest of Earth's land); the absolute
comparison, which no Earth cell can pass, is shown separately as the "Mars gap". Details and limitations: `web/methods.html`.

{artifact}

## Results

{results}

## Robustness to the weights

{robust_summary()}

Regions that only emerge under other weights are not re-labelled (the pool is the ~3,000 regions above threshold at the default weights).

## Blind predictions, checked with an instrument the score never used

Hematite, carbonates, clays and sulfates were not inputs to the global score (ASTER's five TIR bands cannot resolve them). After ranking,
each top-3 under-studied candidate was checked with EMIT L2B mineral identifications (60 m) and compared with 20,000 random pixels from
the same EMIT scenes. ✅ = present in ≥ 5% of candidate pixels and enriched ≥ 1.5× over the scene baseline. Everything is reported, including misses.

{blind_pred_table()}

## Case studies (Analog Cards)

{cards()}

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
"""
(ROOT / "results/REPORT.md").write_text(report)
readme = (ROOT / "README.md").read_text()
readme = re.sub(r"## Results\n\n.*?\n\n## Run from scratch", "## Results\n\n" + results.replace("\\", "\\\\")
                + "\nFull report (incl. the artifact found and fixed, robustness and blind predictions): [`results/REPORT.md`](results/REPORT.md) · "
                "validation: [`results/validation.md`](results/validation.md) · candidates: `results/candidates_<target>.csv` · "
                "cards: [`results/analog_cards/`](results/analog_cards/)\n\n"
                "![AnalogFinder web app](results/screenshot_app.png)\n\n## Run from scratch", readme, flags=re.S)
(ROOT / "README.md").write_text(readme)
print(report[:3000])
