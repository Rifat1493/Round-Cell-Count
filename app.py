"""Test Cellpose + round/glow filter on first MeOH and first Noco images."""

import csv
import os

os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

from pathlib import Path

import tifffile

from src.filter_cells import filter_glowing_round
from src.predict import load_model, run_cellpose

# --- old batch mode (commented out) ---
from src.batch_csv import run_batch

CSV_PATH = Path("data/test.csv")
IMAGE_DIR = Path("data/PP1912C Wash 3.5.23")
OUT_CSV = Path("outputs/result.csv")
DIAMETER = 81.0

if __name__ == "__main__":
    if not CSV_PATH.exists():
        raise FileNotFoundError(CSV_PATH.resolve())
    if not IMAGE_DIR.exists():
        raise FileNotFoundError(IMAGE_DIR.resolve())
    run_batch(CSV_PATH, IMAGE_DIR, OUT_CSV, diameter=DIAMETER)
# --------------------------------------

# IMAGE_DIR = Path("data/PP1912C Wash 3.5.23")
# OUT_DIR = Path("outputs")
# DIAMETER = 100.0  # measured at full resolution; scaled automatically when resize is set
# RESIZE = 720  # longest side; set None for full-res

# # First MeOH (Count 1) and first Noco (Count 1)
# IMAGES = [
#     IMAGE_DIR / "7533 HT29 24H 20X MeOH MediaKept.tif",
#     IMAGE_DIR / "7539 HT29 24H 20X Noco MediaKept.tif",
# ]

# RND_KW = dict(
#     min_roundness=0.70,
#     min_solidity=0.85,
#     min_glow=2.0,
#     min_area=50.0,
#     max_eccentricity=0.85,
#     ring_px=3,
# )

# if __name__ == "__main__":
#     model = load_model()
#     for path in IMAGES:
#         if not path.exists():
#             raise FileNotFoundError(path.resolve())

#         masks, n_total = run_cellpose(
#             path,
#             out_dir=OUT_DIR,
#             model=model,
#             save_viz=True,
#         )
#         img = tifffile.imread(str(path))
#         _, kept, _ = filter_glowing_round(masks, img, **RND_KW)

#         label = "MeOH" if "MeOH" in path.name else "Noco"
#         print("=" * 60)
#         print(f"{label}: {path.name}")
#         print(f"  Count (all cells): {n_total}")
#         print(f"  Rnd Count (round+glow): {len(kept)}")
#     print("=" * 60)

# --- diameter sweep (full resolution) ---
# IMAGE = Path(
#     r"D:\Office\2026\Round-Cell-Count\data\PP1912C Wash 3.5.23\7548 HT29 24H 20X RA172 5 uM MediaReplaced6H.tif"
# )
# OUT_CSV = Path("outputs/diameter_sweep_7534.csv")

# if __name__ == "__main__":
#     if not IMAGE.exists():
#         raise FileNotFoundError(IMAGE.resolve())

#     OUT_CSV.parent.mkdir(parents=True, exist_ok=True)

#     model = load_model()

#     with open(OUT_CSV, "w", newline="", encoding="utf-8") as f:
#         w = csv.DictWriter(f, fieldnames=["diameter", "count"])
#         w.writeheader()
#         f.flush()

#         for d in range(70, 120,3):
#             _, n = run_cellpose(
#                 IMAGE,
#                 diameter=float(d),
#                 resize=None,
#                 model=model,
#                 save_viz=False,
#                 save_masks=False,
#             )
#             w.writerow({"diameter": d, "count": n})
#             f.flush()
#             print(f"diameter={d}  count={n}")

#     print(f"Wrote {OUT_CSV.resolve()}")
