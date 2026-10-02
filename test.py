"""Count round glowing cells from original image + Cellpose mask."""

from pathlib import Path

import tifffile

from src.filter_cells import filter_glowing_round

# Paths — mask must match the original image size
IMAGE_PATH = Path("data/im3.tif")
MASK_PATH = Path("outputs/im3_masks.tif")

# Tunable thresholds (calibrate on a few labelled cells if needed)
MIN_ROUNDNESS = 0.70
MIN_SOLIDITY = 0.85
MIN_GLOW = 2.0  # ring(outside) - interior; raise if too many false positives
MIN_AREA = 50.0
MAX_ECCENTRICITY = 0.85
RING_PX = 3


if __name__ == "__main__":
    if not IMAGE_PATH.exists():
        raise FileNotFoundError(f"Image not found: {IMAGE_PATH.resolve()}")
    if not MASK_PATH.exists():
        raise FileNotFoundError(f"Mask not found: {MASK_PATH.resolve()}")

    image = tifffile.imread(str(IMAGE_PATH))
    masks = tifffile.imread(str(MASK_PATH))
    if masks.ndim > 2:
        masks = masks.squeeze()

    print(f"image: {IMAGE_PATH}  shape={image.shape}")
    print(f"mask:  {MASK_PATH}  shape={masks.shape}  max_label={int(masks.max())}")

    n_total, kept, feats = filter_glowing_round(
        masks,
        image,
        min_roundness=MIN_ROUNDNESS,
        min_solidity=MIN_SOLIDITY,
        min_glow=MIN_GLOW,
        min_area=MIN_AREA,
        max_eccentricity=MAX_ECCENTRICITY,
        ring_px=RING_PX,
    )

    print(f"total cells: {n_total}")
    print(
        f"round+glow cells: {len(kept)}  "
        f"(roundness>={MIN_ROUNDNESS}, solidity>={MIN_SOLIDITY}, "
        f"glow>={MIN_GLOW}, area>={MIN_AREA}, ecc<={MAX_ECCENTRICITY})"
    )
    if feats:
        g, r, s, k = feats["glow"], feats["roundness"], feats["solidity"], feats["keep"]
        print(
            f"pass roundness: {(r >= MIN_ROUNDNESS).sum()}  "
            f"solidity: {(s >= MIN_SOLIDITY).sum()}  "
            f"glow: {(g >= MIN_GLOW).sum()}"
        )
        if len(kept):
            print(f"glow range kept: {g[k].min():.1f} .. {g[k].max():.1f}")
