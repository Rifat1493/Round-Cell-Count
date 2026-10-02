"""Cellpose-SAM load + predict (optional resize, diameter scaling)."""

import os

os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

from pathlib import Path

import cv2
import numpy as np
import torch
from cellpose import models
from cellpose.io import imread, imsave

from src.visualize import image_resize, save_outlines, save_overlay


def load_model(gpu: bool | None = None):
    if gpu is None:
        gpu = torch.cuda.is_available()
    print(f"torch={torch.__version__}  cuda={gpu}")
    if not gpu:
        print("WARNING: CUDA not available — running on CPU (very slow for cpsam_v2).")
    return models.CellposeModel(gpu=gpu, pretrained_model="cpsam_v2")


def run_cellpose(
    image_path: str | Path,
    out_dir: str | Path = "outputs",
    diameter: float | None = 34,
    resize: int | None = None,
    niter: int = 250,
    flow_threshold: float = 0.4,
    cellprob_threshold: float = 0,
    model=None,
    save_viz: bool = True,
    save_masks: bool = True,
):
    """Run Cellpose. If `resize` is set, longest side is scaled to that size first."""
    image_path = Path(image_path)
    out_dir = Path(out_dir)
    if save_masks or save_viz:
        out_dir.mkdir(parents=True, exist_ok=True)

    img_full = imread(str(image_path))
    if resize is not None:
        img = image_resize(img_full, resize=resize)
        # scale diameter to the resized image so cell size stays consistent
        if diameter is not None:
            scale = max(img.shape[:2]) / max(img_full.shape[:2])
            diameter = float(diameter) * scale
    else:
        img = img_full

    print(
        f"input={img_full.shape}  model_input={img.shape}  "
        f"diameter={diameter}  resize={resize}  file={image_path.name}"
    )

    if model is None:
        model = load_model()
    eval_kw = dict(
        niter=niter,
        flow_threshold=flow_threshold,
        cellprob_threshold=cellprob_threshold,
    )
    if diameter is not None:
        eval_kw["diameter"] = diameter
    masks, _, _ = model.eval(img, **eval_kw)

    h0, w0 = img_full.shape[:2]
    if masks.shape[:2] != (h0, w0):
        masks_full = cv2.resize(
            masks.astype(np.uint16), (w0, h0), interpolation=cv2.INTER_NEAREST
        ).astype(np.uint16)
    else:
        masks_full = masks.astype(np.uint16)

    n_cells = int(masks.max())
    stem = image_path.stem
    suffix = f"_r{resize}" if resize is not None else ""
    if save_masks:
        imsave(str(out_dir / f"{stem}{suffix}_masks.tif"), masks_full)
    if save_viz:
        save_outlines(img, masks, out_dir / f"{stem}{suffix}_outlines.png")
        save_overlay(img, masks, out_dir / f"{stem}{suffix}_overlay.png")

    print(f"cells: {n_cells}")
    return masks_full, n_cells
