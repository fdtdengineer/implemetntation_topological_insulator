from __future__ import annotations

import numpy as np


def bloch_hamiltonian(k: float, v: float, w: float, m: float = 0.0) -> np.ndarray:
    return np.array(
        [
            [m, v + w * np.exp(-1j * k)],
            [v + w * np.exp(1j * k), -m],
        ],
        dtype=complex,
    )


def zak_phases(v: float, w: float, m: float = 0.0, n_k: int = 401) -> np.ndarray:
    k_grid = np.linspace(-np.pi, np.pi, n_k, endpoint=False)
    eigenvectors = []
    for k in k_grid:
        _, vectors = np.linalg.eigh(bloch_hamiltonian(k, v, w, m))
        eigenvectors.append(vectors)
    eigenvectors = np.asarray(eigenvectors)

    phases = np.zeros(2)
    for band in range(2):
        loop = 1.0 + 0.0j
        for n in range(n_k):
            overlap = np.vdot(
                eigenvectors[n, :, band],
                eigenvectors[(n + 1) % n_k, :, band],
            )
            if abs(overlap) < 1e-12:
                raise AssertionError("Adjacent eigenvectors have nearly zero overlap.")
            loop *= overlap / abs(overlap)
        phases[band] = np.angle(loop)
    return phases


def open_chain_hamiltonian(n_cells: int, v: float, w: float) -> np.ndarray:
    matrix = np.zeros((2 * n_cells, 2 * n_cells), dtype=float)
    for cell in range(n_cells):
        a, b = 2 * cell, 2 * cell + 1
        matrix[a, b] = matrix[b, a] = v
        if cell < n_cells - 1:
            next_a = 2 * (cell + 1)
            matrix[b, next_a] = matrix[next_a, b] = w
    return matrix


def assert_close(actual: np.ndarray | float, expected: np.ndarray | float, atol: float = 1e-10) -> None:
    if not np.allclose(actual, expected, atol=atol, rtol=0.0):
        raise AssertionError(f"Expected {expected!r}, got {actual!r}")


def main() -> None:
    sigma_z = np.diag([1.0, -1.0])

    # Pass 1: Hermiticity, analytic eigenvalues, and chiral symmetry.
    for k in np.linspace(-np.pi, np.pi, 31):
        hamiltonian = bloch_hamiltonian(k, v=0.8, w=1.1)
        assert_close(hamiltonian, hamiltonian.conj().T)
        expected = np.sqrt(0.8**2 + 1.1**2 + 2 * 0.8 * 1.1 * np.cos(k))
        assert_close(np.linalg.eigvalsh(hamiltonian), np.array([-expected, expected]))
        assert_close(sigma_z @ hamiltonian + hamiltonian @ sigma_z, np.zeros((2, 2)))

    # A staggered onsite term breaks the anticommutation relation.
    massive = bloch_hamiltonian(0.37, v=0.8, w=1.1, m=0.2)
    if np.allclose(sigma_z @ massive + massive @ sigma_z, 0.0, atol=1e-10):
        raise AssertionError("The m term should break chiral symmetry.")

    # Pass 2: Wilson-loop Zak phases in both gapped phases.
    trivial = zak_phases(v=1.2, w=1.0)
    topological = zak_phases(v=0.8, w=1.0)
    assert_close(np.sin(trivial), np.zeros(2), atol=1e-8)
    assert_close(np.cos(trivial), np.ones(2), atol=1e-8)
    assert_close(np.sin(topological), np.zeros(2), atol=1e-8)
    assert_close(np.cos(topological), -np.ones(2), atol=1e-8)

    # Pass 3: finite-chain bulk-edge check.  The topological chain has two
    # exponentially split edge modes, while the trivial chain does not.
    topological_spectrum = np.linalg.eigvalsh(open_chain_hamiltonian(24, v=0.5, w=1.0))
    trivial_spectrum = np.linalg.eigvalsh(open_chain_hamiltonian(24, v=1.5, w=1.0))
    if np.count_nonzero(np.abs(topological_spectrum) < 1e-6) != 2:
        raise AssertionError("Expected two near-zero edge modes in the topological chain.")
    if np.count_nonzero(np.abs(trivial_spectrum) < 1e-3) != 0:
        raise AssertionError("Unexpected near-zero modes in the trivial chain.")

    print("SSH verification passed: algebra, Zak phase, and finite-chain spectra.")


if __name__ == "__main__":
    main()
