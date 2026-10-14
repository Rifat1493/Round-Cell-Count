"""Cellpose counts vs manual Excel counts: DAPI = total cells, PH3 = round (mitotic) cells."""

from pathlib import Path

import numpy as np
import pandas as pd

from src.dataset import image_pairs
from src.imaging import display_gray, load_rgb, save_overlay
from src.segmentation import load_model, reject_dim_regions, segment

DAPI_COLOR = (0, 255, 255)
PH3_COLOR = (0, 255, 0)
COLUMNS = ["folder", "condition", "file_prefix", "replicate", "name_mismatch",
           "dapi_actual", "dapi_pred", "ph3_actual", "ph3_pred"]


def mae(actual, predicted):
    return float(np.mean(np.abs(np.asarray(predicted) - np.asarray(actual))))


def count_one(model, path, folder, out_dir, color_bgr, min_brightness):
    """Segment one image, reject dim regions, save the overlay; returns (count, rejected)."""
    img = load_rgb(path)
    labels, segmented = segment(model, img)
    count = segmented
    if min_brightness is not None:
        labels, count = reject_dim_regions(display_gray(img), labels, min_brightness)
    save_overlay(out_dir / "segmented" / folder / f"{path.stem}_overlay.png", img, labels, color_bgr)
    return count, segmented - count


def write_csv(table, out_dir):
    csv_path = out_dir / "count_evaluation.csv"
    try:
        table.to_csv(csv_path, index=False)
    except PermissionError:
        csv_path = out_dir / f"count_evaluation_{pd.Timestamp.now():%Y%m%d_%H%M%S}.csv"
        table.to_csv(csv_path, index=False)
        print("count_evaluation.csv is open in another program; saved a new copy instead.")
    return csv_path


def evaluate(root, xlsx_path, out_dir, max_pairs=None, ph3_min_brightness=50):
    """max_pairs stops after that many DAPI/PH3 pairs are counted; None runs all of them.

    ph3_min_brightness: PH3 regions whose mean brightness (0-255) is below this are rejected;
    None disables the cutoff. DAPI is never filtered.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    model = load_model()

    rows = []
    missing = []
    for entry, dapi_path, ph3_path in image_pairs(root, xlsx_path):
        if dapi_path is None:
            missing.append(f"{entry['folder']} | {entry['condition']} | rep {entry['replicate']}")
            continue
        dapi_pred, dapi_rejected = count_one(model, dapi_path, entry["folder"], out_dir,
                                             DAPI_COLOR, None)
        ph3_pred, ph3_rejected = count_one(model, ph3_path, entry["folder"], out_dir,
                                           PH3_COLOR, ph3_min_brightness)
        rows.append({**entry, "dapi_pred": dapi_pred, "ph3_pred": ph3_pred})
        print(f"{entry['folder']:28s} {entry['file_prefix']:18s} rep {entry['replicate']}  "
              f"DAPI {entry['dapi_actual']:4d} vs {dapi_pred:4d}   "
              f"PH3 {entry['ph3_actual']:4d} vs {ph3_pred:4d} (-{ph3_rejected})")
        if max_pairs is not None and len(rows) >= max_pairs:
            print(f"\nStopped after {max_pairs} pairs (MAX_PAIRS).")
            break

    table = pd.DataFrame(rows, columns=COLUMNS)
    print("\n" + "=" * 70)
    print(f"Matched pairs: {len(table)}   name mismatches: {int(table['name_mismatch'].sum())}   "
          f"missing: {len(missing)}")
    for item in missing:
        print(f"  missing: {item}")
    if not table.empty:
        print(f"\nDAPI MAE: {mae(table['dapi_actual'], table['dapi_pred']):.2f}")
        print(f"PH3  MAE: {mae(table['ph3_actual'], table['ph3_pred']):.2f}")

    csv_path = write_csv(table, out_dir)
    print(f"\nPer-image table: {csv_path}")
    print(f"Overlays:        {out_dir / 'segmented'}")
    return table
