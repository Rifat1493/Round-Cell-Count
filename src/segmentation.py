"""Cellpose segmentation, then a per-image brightness cutoff on the segmented regions."""

import os

os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

import numpy as np
import torch
from cellpose import models


def load_model(gpu=None):
    if gpu is None:
        gpu = torch.cuda.is_available()
    print(f"torch={torch.__version__}  cuda={gpu}")
    if not gpu:
        print("WARNING: CUDA not available - Cellpose-SAM on CPU is very slow.")
    return models.CellposeModel(gpu=gpu, pretrained_model="cpsam_v2")


def segment(model, img):
    """Run Cellpose on img and return (labels, count)."""
    masks, flows, styles = model.eval(img, cellprob_threshold=2.0)
    labels = masks.astype(np.int32)
    count = len(np.unique(labels[labels > 0]))
    return labels, count


def region_brightness(signal, labels):
    """Mean signal value of every region; returns (region ids, brightness per id)."""
    sums = np.bincount(labels.ravel(), weights=signal.ravel().astype(np.float64))
    areas = np.bincount(labels.ravel())
    ids = np.nonzero(areas)[0]
    ids = ids[ids > 0]
    return ids, sums[ids] / areas[ids]


def reject_dim_regions(signal, labels, min_brightness):
    """Drop regions whose mean brightness (0-255 pixel scale) is below min_brightness.

    signal is the 2D brightness image the regions were found in. Returns (relabelled labels, kept count).
    """
    ids, brightness = region_brightness(signal, labels)
    kept = ids[brightness >= min_brightness]
    lookup = np.zeros(int(labels.max()) + 1, np.int32)
    lookup[kept] = np.arange(1, len(kept) + 1)
    return lookup[labels], len(kept)
