"""Link the manual Excel counts to the DAPI / PH3 image files in the Immunofluorescence folders."""

import re
from pathlib import Path

import pandas as pd

SHEET = "PlantPartsaTubulinIF"
FOLDER_PREFIX = "PP1910B_aTubulinIF_"
DAPI_COLS = (1, 2, 3)
PH3_COLS = (7, 8, 9)
COMPOUND_CODE = re.compile(r"PP\d{4}")


def read_manual_counts(xlsx_path):
    """One dict per (folder, condition, replicate) with the manual DAPI and PH3 counts."""
    sheet = pd.read_excel(xlsx_path, sheet_name=SHEET, header=None)
    entries = []
    folder = None
    for row in sheet.itertuples(index=False):
        name = row[0]
        if not isinstance(name, str):
            continue
        name = name.strip()
        if name.startswith("Overall"):
            break
        if name.startswith(FOLDER_PREFIX):
            folder = name
            continue
        if folder is None:
            continue
        for replicate, (dapi_col, ph3_col) in enumerate(zip(DAPI_COLS, PH3_COLS), 1):
            dapi, ph3 = row[dapi_col], row[ph3_col]
            if pd.isna(dapi) or pd.isna(ph3):
                continue
            entries.append({
                "folder": folder,
                "condition": name,
                "replicate": replicate,
                "dapi_actual": int(dapi),
                "ph3_actual": int(ph3),
            })
    return entries


def condition_to_prefix(condition):
    """'PP-1910B 150 ug/mL' -> 'PP1910B_150ugmL' (compound, then dose with unit joined)."""
    tokens = condition.replace("-", "").replace("/", "").split()
    if len(tokens) == 1:
        return tokens[0]
    return f"{tokens[0]}_{''.join(tokens[1:])}"


def folder_prefixes(folder_path):
    """All condition prefixes in a folder, taken from '<prefix>_DAPI_<n>.tif' files."""
    prefixes = set()
    for path in Path(folder_path).glob("*_DAPI_*.tif"):
        match = re.match(r"(.+)_DAPI_\d+$", path.stem)
        if match:
            prefixes.add(match.group(1))
    return prefixes


def resolve_prefix(folder_path, condition):
    """Return (file_prefix, name_mismatch). Falls back to the PPxxxx compound code."""
    available = folder_prefixes(folder_path)
    prefix = condition_to_prefix(condition)
    if prefix in available:
        return prefix, False
    code = COMPOUND_CODE.search(prefix)
    if code:
        candidates = sorted(p for p in available if p.startswith(code.group(0)))
        if candidates:
            return candidates[0], True
    return None, False


def image_pairs(root, xlsx_path):
    """Yield (entry, dapi_path, ph3_path); paths are None when the files are missing."""
    root = Path(root)
    for entry in read_manual_counts(xlsx_path):
        folder_path = root / entry["folder"]
        prefix, mismatch = resolve_prefix(folder_path, entry["condition"])
        entry = {**entry, "file_prefix": prefix, "name_mismatch": mismatch}
        if prefix is None:
            yield entry, None, None
            continue
        n = entry["replicate"]
        dapi_path = folder_path / f"{prefix}_DAPI_{n}.tif"
        ph3_path = folder_path / f"{prefix}_PH3_{n}.tif"
        if not dapi_path.exists() or not ph3_path.exists():
            yield entry, None, None
            continue
        yield entry, dapi_path, ph3_path
