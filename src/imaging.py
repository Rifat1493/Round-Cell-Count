"""Image loading and overlay drawing."""

from pathlib import Path

import cv2
import numpy as np
import tifffile


def load_rgb(path):
    """Read a TIFF and drop the alpha layer, returning (H, W, 3) or a 2D grayscale image."""
    img = tifffile.imread(str(path))
    if img.ndim == 3 and img.shape[-1] == 4:
        img = img[..., :3]
    return img


def display_gray(img):
    """Single-layer view for drawing: the brightest colour layer at each pixel."""
    return img.max(axis=-1) if img.ndim == 3 else img


def draw_overlay(img, labels, color_bgr):
    """Numbered contours of each label on a contrast-stretched copy of the image."""
    gray = display_gray(img).astype(np.float32)
    low, high = np.percentile(gray, [1, 99.8])
    stretched = np.clip((gray - low) * 255.0 / max(high - low, 1.0), 0, 255).astype(np.uint8)
    out = cv2.cvtColor(stretched, cv2.COLOR_GRAY2BGR)

    for label in np.unique(labels[labels > 0]):
        mask = (labels == label).astype(np.uint8)
        contours, hierarchy = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
        cv2.drawContours(out, contours, -1, color_bgr, 1)
        ys, xs = np.nonzero(mask)
        cv2.putText(out, str(label), (int(xs.mean()) - 6, int(ys.mean()) + 4),
                    cv2.FONT_HERSHEY_SIMPLEX, 0.35, (255, 255, 255), 1, cv2.LINE_AA)
    return out


def save_overlay(path, img, labels, color_bgr):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    cv2.imwrite(str(path), draw_overlay(img, labels, color_bgr))
