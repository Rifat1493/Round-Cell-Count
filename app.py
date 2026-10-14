"""Count cells in the Immunofluorescence images with Cellpose and compare to the manual Excel counts.

DAPI images give the total cell count, PH3 images the round (mitotic) count.
Run in the GPU environment:  conda activate cellpose  &&  python app.py
"""

from pathlib import Path

from src.evaluation import evaluate

IF_ROOT = Path("data/Immunofluorescence")
XLSX_PATH = IF_ROOT / "Arnica Cell Counting IF.xlsx"
OUT_DIR = Path("outputs/evaluation")
MAX_PAIRS = None  # number of DAPI/PH3 image pairs to count; None = all 54
PH3_MIN_BRIGHTNESS = 50  # reject PH3 regions with mean brightness (0-255) below this; None = off

if __name__ == "__main__":
    if not XLSX_PATH.exists():
        raise FileNotFoundError(XLSX_PATH.resolve())
    evaluate(IF_ROOT, XLSX_PATH, OUT_DIR, max_pairs=MAX_PAIRS, ph3_min_brightness=PH3_MIN_BRIGHTNESS)
