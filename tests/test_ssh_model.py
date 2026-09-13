import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from ssh_model import (
    bloch_hamiltonian,
    open_chain_hamiltonian,
    winding_number,
    zak_phases,
)


def test_bloch_hamiltonian_matches_analytic_dispersion() -> None:
    v, w = 0.8, 1.1
    for k in np.linspace(-np.pi, np.pi, 31):
        hamiltonian = bloch_hamiltonian(k, v=v, w=w)
        expected = np.sqrt(v**2 + w**2 + 2 * v * w * np.cos(k))

        assert np.allclose(hamiltonian, hamiltonian.conj().T)
        assert np.allclose(np.linalg.eigvalsh(hamiltonian), [-expected, expected])


def test_chiral_symmetry_and_its_breaking() -> None:
    sigma_z = np.diag([1.0, -1.0])
    chiral = bloch_hamiltonian(0.37, v=0.8, w=1.1)
    massive = bloch_hamiltonian(0.37, v=0.8, w=1.1, mass=0.2)

    assert np.allclose(sigma_z @ chiral + chiral @ sigma_z, 0.0)
    assert not np.allclose(sigma_z @ massive + massive @ sigma_z, 0.0)


def test_winding_number_distinguishes_both_phases() -> None:
    assert winding_number(v=1.2, w=1.0) == 0
    assert winding_number(v=0.8, w=1.0) == 1

    with pytest.raises(ValueError, match="gap closes"):
        winding_number(v=1.0, w=1.0)


def test_zak_phase_is_quantized_in_chiral_model() -> None:
    trivial = zak_phases(v=1.2, w=1.0)
    topological = zak_phases(v=0.8, w=1.0)

    assert np.allclose(np.sin(trivial), 0.0, atol=1e-8)
    assert np.allclose(np.cos(trivial), 1.0, atol=1e-8)
    assert np.allclose(np.sin(topological), 0.0, atol=1e-8)
    assert np.allclose(np.cos(topological), -1.0, atol=1e-8)


def test_bulk_edge_correspondence_for_finite_chain() -> None:
    topological = np.linalg.eigvalsh(open_chain_hamiltonian(24, v=0.5, w=1.0))
    trivial = np.linalg.eigvalsh(open_chain_hamiltonian(24, v=1.5, w=1.0))

    assert np.count_nonzero(np.abs(topological) < 1e-6) == 2
    assert np.count_nonzero(np.abs(trivial) < 1e-3) == 0
