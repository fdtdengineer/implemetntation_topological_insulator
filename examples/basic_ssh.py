"""Minimal usage example for the SSH model package."""

import numpy as np

from ssh_model import open_chain_hamiltonian, winding_number, zak_phases


def main() -> None:
    v, w = 0.8, 1.0
    phases = zak_phases(v=v, w=w)
    winding = winding_number(v=v, w=w)

    finite_chain = open_chain_hamiltonian(n_cells=24, v=v, w=w)
    spectrum = np.linalg.eigvalsh(finite_chain)
    edge_energies = spectrum[np.argsort(np.abs(spectrum))[:2]]

    print(f"winding number: {winding}")
    print(f"Zak phases / pi: {phases / np.pi}")
    print(f"two eigenvalues closest to zero: {edge_energies}")


if __name__ == "__main__":
    main()
