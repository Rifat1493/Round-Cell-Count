"""Map test.csv conditions to PP1912C Wash images and write result.csv."""

import csv
import re
from pathlib import Path

import tifffile

from src.filter_cells import filter_glowing_round
from src.predict import load_model, run_cellpose

# CSV condition name -> filename drug token when they differ
CONDITION_ALIASES = {
    "PP1912B 150 ugmL": "PP1912C 50 ugmL",
    "PP1913B 150 ugmL": "PP1913C 150 ugmL",
}

RND_KW = dict(
    min_roundness=0.70,
    min_solidity=0.85,
    min_glow=2.0,
    min_area=50.0,
    max_eccentricity=0.85,
    ring_px=3,
)


def _numeric_prefix(path: Path) -> int:
    m = re.match(r"^(\d+)", path.name)
    return int(m.group(1)) if m else 10**9


def list_tifs(image_dir: Path) -> list[Path]:
    return sorted(image_dir.glob("*.tif"), key=_numeric_prefix) + sorted(
        image_dir.glob("*.tiff"), key=_numeric_prefix
    )


def files_for_condition(condition: str, tifs: list[Path]) -> list[Path]:
    """MeOH -> MediaKept x3; MeOH wash -> MediaReplaced6H x3 (sorted by file id)."""
    cond = condition.strip()
    if not cond:
        return []
    is_wash = cond.lower().endswith(" wash")
    base = cond[: -len(" wash")].strip() if is_wash else cond
    base = CONDITION_ALIASES.get(base, base)
    media = "MediaReplaced6H" if is_wash else "MediaKept"
    needle = base.lower()
    hits = [
        p
        for p in tifs
        if media.lower() in p.name.lower() and needle in p.name.lower()
    ]
    hits.sort(key=_numeric_prefix)
    return hits[:3]


def read_test_csv(csv_path: Path) -> tuple[str, list[dict]]:
    with open(csv_path, newline="", encoding="utf-8-sig") as f:
        rows = list(csv.reader(f, delimiter="\t"))
    if not rows:
        raise ValueError(f"Empty CSV: {csv_path}")
    header = rows[0]
    experiment = header[0].strip() or "experiment"
    out = []
    for row in rows[1:]:
        if not row or not str(row[0]).strip():
            continue
        # pad short rows
        while len(row) < len(header):
            row.append("")
        out.append({header[i]: row[i] for i in range(len(header))})
    return experiment, out


def predict_image(path: Path, model, out_dir: Path, diameter: float):
    masks, n_total = run_cellpose(
        path, out_dir=out_dir, diameter=diameter, model=model, save_viz=True
    )
    img = tifffile.imread(str(path))
    _, kept, _ = filter_glowing_round(masks, img, **RND_KW)
    return n_total, len(kept)


def run_batch(
    csv_path: Path,
    image_dir: Path,
    out_csv: Path,
    out_dir: Path = Path("outputs"),
    diameter: float = 34.0,
):
    experiment, records = read_test_csv(csv_path)
    tifs = list_tifs(image_dir)
    if not tifs:
        raise FileNotFoundError(f"No .tif files in {image_dir}")

    model = load_model()
    result_rows = []

    for rec in records:
        cond = str(rec.get(experiment, "")).strip()
        paths = files_for_condition(cond, tifs)
        print(f"\n=== {cond} -> {[p.name for p in paths]}")

        pred_counts = ["", "", ""]
        rnd_counts = ["", "", ""]
        file_paths = ["", "", ""]

        for i, p in enumerate(paths):
            file_paths[i] = str(p.resolve())
            try:
                n_tot, n_rnd = predict_image(p, model, out_dir, diameter)
                pred_counts[i] = n_tot
                rnd_counts[i] = n_rnd
            except Exception as e:
                print(f"ERROR {p.name}: {e}")
                pred_counts[i] = f"ERR:{e}"
                rnd_counts[i] = f"ERR:{e}"

        if len(paths) < 3:
            print(f"WARNING: expected 3 images for '{cond}', got {len(paths)}")

        row = dict(rec)
        row["Count 1 predicted"] = pred_counts[0]
        row["Count 2 predicted"] = pred_counts[1]
        row["Count 3 predicted"] = pred_counts[2]
        row["Rnd Count 1 predicted"] = rnd_counts[0]
        row["Rnd Count 2 predicted"] = rnd_counts[1]
        row["Rnd Count 3 predicted"] = rnd_counts[2]
        row["Count 1 file path"] = file_paths[0]
        row["Count 2 file path"] = file_paths[1]
        row["Count 3 file path"] = file_paths[2]
        result_rows.append(row)

    # Fixed column order: Count N next to Count N predicted
    fieldnames = [experiment]
    for n in (1, 2, 3):
        fieldnames += [f"Count {n}", f"Count {n} predicted"]
    fieldnames += ["Average"]
    for n in (1, 2, 3):
        fieldnames += [f"Rnd Count {n}", f"Rnd Count {n} predicted"]
    fieldnames += ["Count 1 file path", "Count 2 file path", "Count 3 file path"]
    # keep any unexpected original columns
    for r in result_rows:
        for k in r:
            if k not in fieldnames:
                fieldnames.append(k)

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    # Comma-separated + UTF-8 BOM so Excel shows real columns
    with open(out_csv, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter=",", extrasaction="ignore")
        w.writeheader()
        w.writerows(result_rows)
    print(f"\nWrote {out_csv}")
    return out_csv
