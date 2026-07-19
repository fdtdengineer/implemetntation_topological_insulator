from __future__ import annotations

import json
import re
from pathlib import Path

path = Path("sshmodel.ipynb")
notebook = json.loads(path.read_text(encoding="utf-8"))


def get_text(cell):
    source = cell.get("source", [])
    return "".join(source) if isinstance(source, list) else str(source)


def put_text(cell, text):
    cell["source"] = text.splitlines(keepends=True)


H_BLOCH = '''def H_SSH(k, v, w, m=0.0):
    """Bloch Hamiltonian in the (A, B) sublattice basis."""
    return np.array(
        [[m, v + w * np.exp(-1j * k)],
         [v + w * np.exp(1j * k), -m]],
        dtype=complex,
    )
'''

H_OPEN = '''def H_SSH_open(n_cells, v, w):
    """Open SSH chain in the basis (A1, B1, A2, B2, ...)."""
    H = np.zeros((2 * n_cells, 2 * n_cells), dtype=float)
    for cell in range(n_cells):
        a, b = 2 * cell, 2 * cell + 1
        H[a, b] = H[b, a] = v
        if cell < n_cells - 1:
            next_a = 2 * (cell + 1)
            H[b, next_a] = H[next_a, b] = w
    return H
'''

FINITE = '''L = 20
v_values = np.linspace(0.0, 3.0, 61)
w = 1.0
spectra = []
eigenvectors_open = []

for v in v_values:
    values, vectors = np.linalg.eigh(H_SSH_open(L, v, w))
    spectra.append(values)
    eigenvectors_open.append(vectors)

spectra = np.asarray(spectra)
eigenvectors_open = np.asarray(eigenvectors_open)
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

    probability = np.abs(eigenvectors_open[vi, :, n]) ** 2
    plt.subplot(133)
    plt.cla()
    plt.bar(np.arange(0, 2 * L, 2), probability[0::2], 0.9, label="A")
    plt.bar(np.arange(1, 2 * L, 2), probability[1::2], 0.9, label="B")
    plt.ylim(0.0, 1.0)
    plt.yticks(np.linspace(0.0, 1.0, 5))
    plt.ylabel(r"$|\\psi|^2$", fontsize=fs)
    plt.xlabel("site index", fontsize=fs)
    plt.legend(loc="upper right")
'''

GEN_DATA = '''def gen_data(H, v, w, m=0.0, N=101):
    """Diagonalize a gapped Hermitian Bloch Hamiltonian on [-pi, pi)."""
    k_list = np.linspace(-np.pi, np.pi, N, endpoint=False)
    eigenvalues = []
    eigenvectors = []
    for k in k_list:
        values, vectors = np.linalg.eigh(H(k, v, w, m=m))
        eigenvalues.append(values)
        eigenvectors.append(vectors)
    return k_list, np.asarray(eigenvalues), np.asarray(eigenvectors)
'''

RUN = '''w = 1.0
v = 0.9
m = 0.0
N = 101
k_list, eigenvalues, eigenvectors = gen_data(H_SSH, v, w, m=m, N=N)
# bnd_plot(k_list, eigenvalues, bandtype="Re")
'''

BAND_PLOT = '''def bnd_plot(k_list, eigenvalues, bandtype="Re", fs=18):
    plt.figure(figsize=(5, 4))
    values = eigenvalues.real if bandtype == "Re" else eigenvalues.imag
    for band in range(values.shape[1]):
        plt.plot(k_list, values[:, band])
    plt.xlabel(r"$k_xa$", fontsize=fs)
    plt.ylabel(f"{bandtype}$(E)$", fontsize=fs)
    plt.xticks(
        np.linspace(-np.pi, np.pi, 5),
        [r"$-\\pi$", r"$-\\pi/2$", r"$0$", r"$\\pi/2$", r"$\\pi$"],
        fontsize=fs,
    )
    plt.tick_params(labelsize=0.8 * fs)
    plt.show()
'''

WILSON = '''def zak_phase_from_eigenvectors(eigenvectors, atol=1e-12):
    """Return one gauge-invariant Zak phase per isolated band."""
    n_k, _, n_band = eigenvectors.shape
    phases = np.zeros(n_band)
    for band in range(n_band):
        wilson_loop = 1.0 + 0.0j
        for n in range(n_k):
            overlap = np.vdot(
                eigenvectors[n, :, band],
                eigenvectors[(n + 1) % n_k, :, band],
            )
            if abs(overlap) < atol:
                raise ValueError("Adjacent eigenvectors have nearly zero overlap.")
            wilson_loop *= overlap / abs(overlap)
        phases[band] = np.angle(wilson_loop)
    return phases

zak_phase = zak_phase_from_eigenvectors(eigenvectors)
print("Zak phase / pi:", zak_phase / np.pi)
'''

new_cells = []
for cell in notebook.get("cells", []):
    text = get_text(cell)

    if cell.get("cell_type") == "code" and any(
        marker in text
        for marker in ("npr_psiL[:,:,0]", "npr_psiL[:, :, 0]", "phase_pi = np.angle(phase)")
    ):
        continue

    if cell.get("cell_type") == "markdown":
        text = text.replace(
            "うろ覚えでかなりいい加減なことを書いているので、詳しくはその手のちゃんとしたテキストを参照すること。",
            "以下では、標準的な単位胞と実数結合 $v,w>0$ を基準に、対称性とトポロジカル不変量を整理する。",
        )
        text = re.sub(
            r"上の実装では、追加でChiral対称性を壊す項.*?後ほど確認する。",
            "上の実装では、カイラル対称性を破る staggered onsite term $m\\sigma_z$ も扱える。$m\\neq0$ では Zak phase は一般に $0$ または $\\pi$ に量子化されず、零エネルギー端状態もカイラル対称性によって保護されなくなる。",
            text,
        )
        text = text.replace(
            "一般の符号を許す場合は $\\lvert w\\rvert$ と $\\lvert v\\rvert$ の大小に加えて、曲線の向きによって $\\nu$ の符号も変わり得る。",
            "この Fourier 規約では、実数結合について $|w|>|v|$ なら $\\nu=1$、$|w|<|v|$ なら $\\nu=0$ である。逆向きの Fourier 規約を採用すると winding number の符号は反転するが、相の分類は変わらない。",
        )
        text = text.replace(
            "上のコードを適当に弄ると、非エルミートな場合の計算も可能。各自試してみてください。",
            "",
        )
        put_text(cell, text)

    elif cell.get("cell_type") == "code":
        if "def H_SSH_reals" in text or "def H_SSH_open" in text:
            text = H_OPEN
        elif "L=20" in text and "def enpsi" in text:
            text = FINITE
        elif "def H_SSH(" in text:
            text = H_BLOCH
        elif "def gen_data(" in text:
            text = GEN_DATA
        elif "k_list, eigenvalues" in text and "gen_data(" in text:
            text = RUN
        elif "def bnd_plot(" in text:
            text = BAND_PLOT
        elif "wilson_loop = np.ones" in text or "def zak_phase_from_eigenvectors" in text:
            text = WILSON
        put_text(cell, text)
        cell["outputs"] = []
        cell["execution_count"] = None

    new_cells.append(cell)

notebook["cells"] = new_cells
path.write_text(json.dumps(notebook, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")
