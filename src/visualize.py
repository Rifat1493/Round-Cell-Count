"""Image helpers for Cellpose outputs."""

import cv2
import numpy as np


def normalize99(img: np.ndarray) -> np.ndarray:
    x = img.astype(np.float32)
    lo, hi = np.percentile(x, 1), np.percentile(x, 99)
    return (x - lo) / (1e-10 + hi - lo)


def image_resize(img: np.ndarray, resize: int = 720) -> np.ndarray:
    """Shrink so the longest side equals `resize` (no upscale)."""
    ny, nx = img.shape[:2]
    if max(ny, nx) <= resize:
        return img
    if ny > nx:
        nx, ny = int(nx / ny * resize), resize
    else:
        ny, nx = int(ny / nx * resize), resize
    return cv2.resize(img, (nx, ny))


def save_outlines(img: np.ndarray, masks: np.ndarray, path) -> None:
    gray = normalize99(img.mean(axis=-1) if img.ndim == 3 else img)
    base = cv2.cvtColor(np.clip(gray * 255, 0, 255).astype(np.uint8), cv2.COLOR_GRAY2BGR)
    contours, _ = cv2.findContours(
        masks.astype(np.int32), cv2.RETR_FLOODFILL, cv2.CHAIN_APPROX_SIMPLE
    )
    keep = [c for c in contours if len(c.squeeze().shape) == 2 and len(c) > 4]
    cv2.drawContours(base, keep, -1, (0, 0, 255), 1)
    cv2.imwrite(str(path), base)


def save_overlay(img: np.ndarray, masks: np.ndarray, path) -> None:
    gray = normalize99(img.mean(axis=-1) if img.ndim > 2 else img)
    hsv = np.zeros((*gray.shape, 3), np.uint8)
    hsv[:, :, 2] = np.clip(gray * 1.5 * 255, 0, 255).astype(np.uint8)
    for n in range(int(masks.max())):
        ipix = (masks == n + 1).nonzero()
        hsv[ipix[0], ipix[1], 0] = np.random.randint(0, 180)
        hsv[ipix[0], ipix[1], 1] = 255
    cv2.imwrite(str(path), cv2.cvtColor(hsv, cv2.COLOR_HSV2BGR))
