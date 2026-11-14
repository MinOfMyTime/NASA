"""Sanity-check the Mars thermal model against Spanovich et al. (2006): Opportunity max surface T = 295 K at Ls=339."""
import json, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from mars_thermal_model import diurnal_cycle  # noqa: E402
fp = json.load(open(Path(__file__).resolve().parents[1] / "data/processed/mars_fingerprints.json"))["meridiani"]
ti, alb = fp["features"]["thermal_inertia"]["median"], fp["features"]["albedo"]["median"]
tmin, tmax = diurnal_cycle(ti, alb, -1.95, 339)
print(f"Meridiani TI={ti:.0f} albedo={alb:.3f}: model Tmax={tmax:.1f} K Tmin={tmin:.1f} K (observed Tmax 295 K)")
