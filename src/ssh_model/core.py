"""Core numerical routines for the Su-Schrieffer-Heeger model."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

ComplexArray = NDArray[np.complex128]
FloatArray = NDArray[np.float64]


def bloch_hamiltonian(k: float, v: float, w: float, mass: float = 0.0) -> ComplexArray:
    """Return the 2x2 SSH Bloch Hamiltonian in the (A, B) basis."""
    hopping = v + w * np.exp(1j * k)
    return np.array(
        [[mass, np.conjugate(hopping)], [hopping, -mass]],
        dtype=np.complex128,
    )


def band_structure(
    v: float,
    w: float,
    mass: float = 0.0,
    n_k: int = 401,
) -> tuple[FloatArray, FloatArray, ComplexArray]:
    """Compute eigenvalues and eigenvectors on a uniform Brillouin-zone grid."""
    if n_k < 2:
        raise ValueError("n_k must be at least 2.")

    k_grid = np.linspace(-np.pi, np.pi, n_k, endpoint=False)
    eigenvalues = np.empty((n_k, 2), dtype=float)
    eigenvectors = np.empty((n_k, 2, 2), dtype=np.complex128)

    for index, k in enumerate(k_grid):
        values, vectors = np.linalg.eigh(bloch_hamiltonian(k, v, w, mass))
        eigenvalues[index] = values
        eigenvectors[index] = vectors

    return k_grid, eigenvalues, eigenvectors


def zak_phase_from_eigenvectors(
    eigenvectors: ComplexArray,
    atol: float = 1e-12,
) -> FloatArray:
    """Return one gauge-invariant Wilson-loop Zak phase per isolated band."""
    if eigenvectors.ndim != 3 or eigenvectors.shape[1:] != (2, 2):
        raise ValueError("eigenvectors must have shape (n_k, 2, 2).")

    n_k = eigenvectors.shape[0]
    phases = np.zeros(2, dtype=float)

    for band in range(2):
        wilson_loop = 1.0 + 0.0j
        for index in range(n_k):
            state = eigenvectors[index, :, band]
            next_state = eigenvectors[(index + 1) % n_k, :, band]
            overlap = np.vdot(state, next_state)
            if abs(overlap) < atol:
                raise ValueError("Adjacent eigenvectors have nearly zero overlap.")
            wilson_loop *= overlap / abs(overlap)
        phases[band] = np.angle(wilson_loop)

    return phases


def zak_phases(
    v: float,
    w: float,
    mass: float = 0.0,
    n_k: int = 401,
) -> FloatArray:
    """Compute the Wilson-loop Zak phase for both SSH bands."""
    _, _, eigenvectors = band_structure(v=v, w=w, mass=mass, n_k=n_k)
    return zak_phase_from_eigenvectors(eigenvectors)


def winding_number(v: float, w: float, n_k: int = 4097, atol: float = 1e-12) -> int:
    """Return the winding number of h(k)=v+w exp(ik) for a gapped SSH model."""
    if n_k < 3:
        raise ValueError("n_k must be at least 3.")

    k_grid = np.linspace(-np.pi, np.pi, n_k)
    h = v + w * np.exp(1j * k_grid)
    if np.min(np.abs(h)) < atol:
        raise ValueError("The bulk gap closes; the winding number is undefined.")

    phase = np.unwrap(np.angle(h))
    return int(np.rint((phase[-1] - phase[0]) / (2 * np.pi)))


def open_chain_hamiltonian(n_cells: int, v: float, w: float) -> FloatArray:
    """Return the real-space Hamiltonian for an open SSH chain."""
    if n_cells < 1:
        raise ValueError("n_cells must be positive.")

    matrix = np.zeros((2 * n_cells, 2 * n_cells), dtype=float)
    for cell in range(n_cells):
        a_site = 2 * cell
        b_site = a_site + 1
        matrix[a_site, b_site] = matrix[b_site, a_site] = v

        if cell < n_cells - 1:
            next_a_site = 2 * (cell + 1)
            matrix[b_site, next_a_site] = matrix[next_a_site, b_site] = w

    return matrix
