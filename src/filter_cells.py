"""Filter Cellpose masks for round, glowing cells (OpenCV + NumPy only)."""

import cv2
import numpy as np


def _to_gray(img: np.ndarray) -> np.ndarray:
    if img.ndim == 3:
        return img.astype(np.float32).mean(axis=-1)
    return img.astype(np.float32)


def _features_for_label(masks: np.ndarray, gray: np.ndarray, label: int, ring_px: int):
    ys, xs = np.where(masks == label)
    if len(ys) == 0:
        return None

    h, w = masks.shape
    y0, y1 = max(0, ys.min() - ring_px), min(h, ys.max() + ring_px + 1)
    x0, x1 = max(0, xs.min() - ring_px), min(w, xs.max() + ring_px + 1)
    crop_m = masks[y0:y1, x0:x1]
    crop_g = gray[y0:y1, x0:x1]
    inner = (crop_m == label).astype(np.uint8)

    contours, _ = cv2.findContours(inner, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    if not contours:
        return None
    c = max(contours, key=cv2.contourArea)
    area = float(cv2.contourArea(c))
    if area <= 0:
        return None

    hull = cv2.convexHull(c)
    hull_area = float(cv2.contourArea(hull)) or 1.0
    solidity = area / hull_area

    (_, _), radius = cv2.minEnclosingCircle(c)
    roundness = area / (np.pi * radius * radius) if radius > 0 else 0.0

    if len(c) >= 5:
        (_, _), (ma, mb), _ = cv2.fitEllipse(c)
        major, minor = max(ma, mb), min(ma, mb)
        eccentricity = float(np.sqrt(max(0.0, 1.0 - (minor / major) ** 2))) if major > 0 else 1.0
    else:
        eccentricity = 1.0

    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (2 * ring_px + 1, 2 * ring_px + 1))
    dil = cv2.dilate(inner, kernel)
    ring = (dil > 0) & (inner == 0) & (crop_m == 0)
    if ring.any() and inner.any():
        glow = float(crop_g[ring].mean() - crop_g[inner.astype(bool)].mean())
    else:
        glow = 0.0

    return {
        "label": label,
        "area": area,
        "roundness": float(roundness),
        "solidity": float(solidity),
        "eccentricity": eccentricity,
        "glow": glow,
    }


def filter_glowing_round(
    masks: np.ndarray,
    image: np.ndarray,
    min_roundness: float = 0.70,
    min_solidity: float = 0.85,
    min_glow: float = 5.0,
    min_area: float = 50.0,
    max_eccentricity: float = 0.85,
    ring_px: int = 3,
):
    """Return (n_total, kept_labels, feature_dict) for round glowing cells."""
    if masks.ndim > 2:
        masks = masks.squeeze()
    gray = _to_gray(image)
    if gray.shape[:2] != masks.shape[:2]:
        raise ValueError(f"image {gray.shape} and masks {masks.shape} size mismatch")

    n_total = int(masks.max())
    if n_total == 0:
        return 0, [], {}

    rows = []
    for lab in range(1, n_total + 1):
        feat = _features_for_label(masks, gray, lab, ring_px)
        if feat is not None:
            rows.append(feat)

    if not rows:
        return n_total, [], {}

    labels = np.array([r["label"] for r in rows])
    area = np.array([r["area"] for r in rows], dtype=np.float32)
    roundness = np.array([r["roundness"] for r in rows], dtype=np.float32)
    solidity = np.array([r["solidity"] for r in rows], dtype=np.float32)
    eccentricity = np.array([r["eccentricity"] for r in rows], dtype=np.float32)
    glow = np.array([r["glow"] for r in rows], dtype=np.float32)

    keep = (
        (area >= min_area)
        & (roundness >= min_roundness)
        & (solidity >= min_solidity)
        & (eccentricity <= max_eccentricity)
        & (glow >= min_glow)
    )
    kept = labels[keep].tolist()
    features = {
        "label": labels,
        "area": area,
        "roundness": roundness,
        "solidity": solidity,
        "eccentricity": eccentricity,
        "glow": glow,
        "keep": keep,
    }
    return n_total, kept, features
