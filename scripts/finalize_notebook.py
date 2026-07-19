from __future__ import annotations

import json
from pathlib import Path

NOTEBOOK_PATH = Path("sshmodel.ipynb")

H_SSH_CODE = r'''def H_SSH(k, v, w, m=0.0):
    """Bloch Hamiltonian in the (A, B) sublattice basis."""
    return np.array(
        [
            [m, v + w * np.exp(-1j * k)],
            [v + w * np.exp(1j * k), -m],
        ],
        dtype=complex,
    )
'''

OPEN_CHAIN_CODE = r'''def H_SSH_open(n_cells, v, w):
    """Open SSH chain with basis (A1, B1, A2, B2, ...)."""
    H = np.zeros((2 * n_cells, 2 * n_cells), dtype=float)

    for cell in range(n_cells):
        a = 2 * cell
        b = a + 1
        H[a, b] = H[b, a] = v

        if cell < n_cells - 1:
            next_a = 2 * (cell + 1)
            H[b, next_a] = H[next_a, b] = w

    return H
'''

FINITE_SWEEP_CODE = r'''L = 20
v_values = np.linspace(0.0, 3.0, 61)
w = 1.0
spectra = []
eigenvectors = []

for v in v_values:
    eigenvalues, vectors = np.linalg.eigh(H_SSH_open(L, v, w))
    spectra.append(eigenvalues)
    eigenvectors.append(vectors)

spectra = np.asarray(spectra)
eigenvectors = np.asarray(eigenvectors)

fs = 15
plt.figure(figsize=(15, 6))

@ipywidgets.interact(vi=(0, len(v_values) - 1), n=(0, 2 * L - 1))
def plot_finite_chain(vi=20, n=L):
    plt.subplot(131)
    plt.cla()
    plt.plot(np.arange(2 * L), spectra[vi], "ko")
    plt.xlabel("eigenstate index", fontsize=fs)
    plt.ylabel(r"$E$", fontsize=fs)
    plt.ylim(-3.0, 3.0)

    plt.subplot(132)
    plt.cla()
    plt.plot(v_values, spectra, "k-", linewidth=0.8)
    plt.plot(v_values[vi], spectra[vi, n], "ro", markersize=8)
    plt.axvline(w, linestyle="--", linewidth=1)
    plt.xlabel(r"$v$", fontsize=fs)
    plt.ylabel(r"$E$", fontsize=fs)
    plt.ylim(-3.0, 3.0)

    probability = np.abs(eigenvectors[vi, :, n]) ** 2
    plt.subplot(133)
    plt.cla()
    plt.bar(np.arange(0, 2 * L, 2), probability[0::2], 0.9, label="A")
    plt.bar(np.arange(1, 2 * L, 2), probability[1::2], 0.9, label="B")
    plt.ylim(0.0, 1.0)
    plt.yticks(np.linspace(0.0, 1.0, 5))
    plt.ylabel(r"$|\psi|^2$", fontsize=fs)
    plt.xlabel("site index", fontsize=fs)
    plt.legend(loc="upper right")
'''

GEN_DATA_CODE = r'''def gen_data(H, v, w, m=0.0, N=101):
    """Diagonalize a gapped Hermitian Bloch Hamiltonian on [-pi, pi)."""
    k_list = np.linspace(-np.pi, np.pi, N, endpoint=False)
    eigenvalues = []
    eigenvectors = []

    for k in k_list:
        E, U = np.linalg.eigh(H(k, v, w, m=m))
        eigenvalues.append(E)
        eigenvectors.append(U)

    return k_list, np.asarray(eigenvalues), np.asarray(eigenvectors)
'''

RUN_CODE = r'''w = 1.0
v = 0.9
m = 0.0
N = 101

k_list, eigenvalues, eigenvectors = gen_data(H_SSH, v, w, m=m, N=N)
# bnd_plot(k_list, eigenvalues, bandtype="Re")
'''

BAND_PLOT_CODE = r'''def bnd_plot(k_list, eigenvalues, bandtype="Re", fs=18):
    plt.figure(figsize=(5, 4))
    values = eigenvalues.real if bandtype == "Re" else eigenvalues.imag

    for band in range(values.shape[1]):
        plt.plot(k_list, values[:, band])

    plt.xlabel(r"$k_xa$", fontsize=fs)
    plt.ylabel(f"{bandtype}$(E)$", fontsize=fs)
    plt.xticks(
        np.linspace(-np.pi, np.pi, 5),
        [r"$-\pi$", r"$-\pi/2$", r"$0$", r"$\pi/2$", r"$\pi$"],
        fontsize=fs,
    )
    plt.tick_params(labelsize=0.8 * fs)
    plt.show()
'''

WILSON_CODE = r'''def zak_phase_from_eigenvectors(eigenvectors, atol=1e-12):
    """Return one gauge-invariant Zak phase per isolated band."""
    n_k, _, n_band = eigenvectors.shape
    phases = np.zeros(n_band)

    for band in range(n_band):
        wilson_loop = 1.0 + 0.0j
        for n in range(n_k):
            u_n = eigenvectors[n, :, band]
            u_next = eigenvectors[(n + 1) % n_k, :, band]
            overlap = np.vdot(u_n, u_next)
            if abs(overlap) < atol:
                raise ValueError("Adjacent eigenvectors have nearly zero overlap.")
            wilson_loop *= overlap / abs(overlap)

        phases[band] = np.angle(wilson_loop)

    return phases

zak_phase = zak_phase_from_eigenvectors(eigenvectors)
print("Zak phase / pi:", zak_phase / np.pi)
'''


def text_of(cell: dict) -> str:
    source = cell.get("source", [])
    return "".join(source) if isinstance(source, list) else str(source)


def set_text(cell: dict, text: str) -> None:
    cell["source"] = text.splitlines(keepends=True)


def should_remove(cell: dict, text: str) -> bool:
    if cell.get("cell_type") == "code":
        stale_markers = (
            "npr_psiL[:,:,0]",
            "npr_psiL[:, :, 0]",
            "phase_pi = np.angle(phase)",
        )
        return any(marker in text for marker in stale_markers)

    if cell.get("cell_type") == "markdown":
        return "上のコードを適当に弄ると、非エルミートな場合の計算も可能" in text

    return False


def rewrite_markdown(text: str) -> str:
    text = text.replace(
        "うろ覚えでかなりいい加減なことを書いているので、詳しくはその手のちゃんとしたテキストを参照すること。\n",
        "以下では、標準的な単位胞と実数結合 $v,w>0$ を基準に対称性とトポロジカル不変量を整理する。\n",
    )
    text = text.replace(
        "上の実装では、追加でChiral対称性を壊す項 $m$ を加えている。$m\\neq 0$ のときにトポロジカル性が失われることを後ほど確認する。",
        "上の実装では、カイラル対称性を破る staggered onsite term $m\\sigma_z$ も扱える。$m\\neq0$ では Zak phase は一般に $0$ または $\\pi$ に量子化されず、SSH の零エネルギー端状態もカイラル対称性によって保護されなくなる。",
    )
    return text


def rewrite_code(text: str) -> str:
    if "def H_SSH_reals" in text or "def H_SSH_open" in text:
        return OPEN_CHAIN_CODE
    if "L=20" in text and "def enpsi" in text:
        return FINITE_SWEEP_CODE
    if "def H_SSH(" in text:
        return H_SSH_CODE
    if "def gen_data(" in text:
        return GEN_DATA_CODE
    if "k_list, eigenvalues" in text and "gen_data(" in text:
        return RUN_CODE
    if "def bnd_plot(" in text:
        return BAND_PLOT_CODE
    if "wilson_loop = np.ones" in text or "def zak_phase_from_eigenvectors" in text:
        return WILSON_CODE
    return text


def validate(notebook: dict) -> None:
    markdown = "\n".join(
        text_of(cell) for cell in notebook["cells"] if cell.get("cell_type") == "markdown"
    )
    code = "\n".join(
        text_of(cell) for cell in notebook["cells"] if cell.get("cell_type") == "code"
    )

    required_markdown = (
        r"\{\sigma_z,H(k)\}=0",
        r"\frac{1}{2}\int_{-\pi}^{\pi}",
        r"\gamma_\pm=\operatorname{Arg}W_\pm",
        r"H_{1,2N}=H_{2N,1}=w",
        r"m\neq0",
    )
    required_code = (
        "np.linalg.eigh",
        "overlap / abs(overlap)",
        "(n + 1) % n_k",
        "probability = np.abs(eigenvectors[vi, :, n]) ** 2",
    )
    forbidden = (
        "phase_pi = np.angle(phase)",
        "npr_psiL[:,:,0]",
        "npr_psiL[:, :, 0]",
        "wilson_loop = np.log(wilson_loop) / (2*np.pi*1.j)",
        "Willson loop",
        "Hermite の場合",
        "plt.yticks(np.linspace(-1,1,5))",
        "上のコードを適当に弄ると、非エルミートな場合の計算も可能",
    )

    for expression in required_markdown:
        if expression not in markdown:
            raise ValueError(f"Required Markdown expression is missing: {expression}")
    for expression in required_code:
        if expression not in code:
            raise ValueError(f"Required code expression is missing: {expression}")
    for expression in forbidden:
        if expression in markdown or expression in code:
            raise ValueError(f"Stale or incorrect content remains: {expression}")


with NOTEBOOK_PATH.open("r", encoding="utf-8") as file:
    notebook = json.load(file)

new_cells = []
for cell in notebook.get("cells", []):
    text = text_of(cell)
    if should_remove(cell, text):
        continue

    if cell.get("cell_type") == "markdown":
        set_text(cell, rewrite_markdown(text))
    elif cell.get("cell_type") == "code":
        set_text(cell, rewrite_code(text))
        cell["outputs"] = []
        cell["execution_count"] = None

    new_cells.append(cell)

notebook["cells"] = new_cells
validate(notebook)

with NOTEBOOK_PATH.open("w", encoding="utf-8") as file:
    json.dump(notebook, file, ensure_ascii=False, indent=1)
    file.write("\n")
