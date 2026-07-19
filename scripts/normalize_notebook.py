from __future__ import annotations

import json
import re
from pathlib import Path

NOTEBOOK_PATH = Path("sshmodel.ipynb")

TITLE = "# 1次元トポロジカル絶縁体：SSH模型の実装"
CITATION = (
    "## Reference\n\n"
    "J. K. Asbóth, L. Oroszlány, and A. Pályi, "
    "*A Short Course on Topological Insulators: Band-Structure Topology and Edge States in One and Two Dimensions*, "
    "Lecture Notes in Physics **919** (Springer, 2016). "
    "[doi:10.1007/978-3-319-25607-8](https://doi.org/10.1007/978-3-319-25607-8); "
    "[arXiv:1509.02295](https://arxiv.org/abs/1509.02295).\n\n"
)

PLANE_WAVE_SECTION = r"""### 平面波基底による表現

周期境界条件の下では、各サブ格子 \(\alpha\in\{A,B\}\) に対して

$$
\lvert k,\alpha\rangle
=\frac{1}{\sqrt{N}}\sum_{l=1}^{N}e^{ikl}\lvert l,\alpha\rangle,
\qquad
k=\frac{2\pi n}{N},\quad n=0,1,\ldots,N-1
$$

と定義する。基底を \((\lvert k,A\rangle,\lvert k,B\rangle)\) の順に取ると、Bloch Hamiltonian は

$$
H(k)=
\begin{bmatrix}
0 & v+w e^{-ik}\\
v+w e^{ik} & 0
\end{bmatrix}
=
\begin{bmatrix}
0 & h(k)^*\\
h(k) & 0
\end{bmatrix},
\qquad
h(k)=v+w e^{ik}
$$

となる。したがって、バンドギャップは \(h(k)=0\) のときに閉じる。標準的な \(v,w>0\) の場合、転移点は \(v=w\)、波数は \(k=\pi\) である。
"""

EIGENVECTOR_SECTION = r"""### 固有値と固有ベクトル

\(h(k)=\lvert h(k)\rvert e^{i\phi(k)}\) と書くと、固有値と規格化固有ベクトルは

$$
E_\pm(k)=\pm\lvert h(k)\rvert
=\pm\sqrt{v^2+w^2+2vw\cos k},
$$

$$
\lvert u_\pm(k)\rangle
=\frac{1}{\sqrt{2}}
\begin{bmatrix}
1\\
\pm e^{i\phi(k)}
\end{bmatrix},
\qquad
\phi(k)=\operatorname{Arg}\!\left(v+w e^{ik}\right)
=\operatorname{atan2}\!\left(w\sin k,\,v+w\cos k\right).
$$

単純な \(\arctan(y/x)\) では象限と分岐を正しく追えないため、数値計算では \(\operatorname{atan2}\) または複素偏角を用いる。ギャップ閉鎖点では \(h(k)=0\) となり、\(\phi(k)\) と上の固有ベクトル表示は定義できない。
"""

CHIRAL_SECTION = r"""### Chiral symmetry

Bloch Hamiltonian は

$$
H(k)=d_x(k)\sigma_x+d_y(k)\sigma_y,
\qquad
d_x(k)=v+w\cos k,\quad d_y(k)=w\sin k,
$$

と書け、\(\sigma_z\) 成分を持たない。このため

$$
\{\sigma_z,H(k)\}=0
\qquad\Longleftrightarrow\qquad
\sigma_zH(k)\sigma_z^{-1}=-H(k)
$$

が成り立つ。これがカイラル対称性である。\(\sigma_z=\operatorname{diag}(1,-1)\) は A/B サブ格子を交換する操作ではなく、両サブ格子に相対的な符号を与える操作である。A/B の交換は \(\sigma_x\) に対応する。
"""

WINDING_SECTION = r"""### Winding number

ギャップが開いている \(h(k)\neq0\) の場合、\(k\) が Brillouin zone を一周するときの複素数 \(h(k)=\lvert h(k)\rvert e^{i\phi(k)}\) の巻き付き数を

$$
\nu
=\frac{1}{2\pi i}\int_{-\pi}^{\pi}
\frac{\mathrm{d}}{\mathrm{d}k}\log h(k)\,\mathrm{d}k
=\frac{1}{2\pi}\int_{-\pi}^{\pi}
\frac{\partial\phi(k)}{\partial k}\,\mathrm{d}k
$$

で定義する。\(\nu\) は、\(h(k)\) の閉曲線が複素平面の原点を向き付きで何周するかを表す整数である。

標準的な実数結合 \(v,w>0\) では

$$
\nu=
\begin{cases}
0, & v>w \quad\text{(trivial)},\\
1, & v<w \quad\text{(topological)},
\end{cases}
$$

であり、\(v=w\) ではギャップが閉じるため winding number は定義できない。一般の符号を許す場合は \(\lvert w\rvert\) と \(\lvert v\rvert\) の大小に加えて、曲線の向きによって \(\nu\) の符号も変わり得る。
"""

ZAK_SECTION = r"""### Zak phase

ここでは Berry connection を

$$
A_\pm(k)
=-i\left\langle u_\pm(k)\middle|\frac{\partial}{\partial k}\middle|u_\pm(k)\right\rangle
$$

と定義し、Zak phase を

$$
\gamma_\pm=\int_{-\pi}^{\pi}A_\pm(k)\,\mathrm{d}k
\qquad (\mathrm{mod}\ 2\pi)
$$

とする。上の固有ベクトルでは

$$
\left\langle u_\pm(k)\middle|\frac{\partial}{\partial\phi}\middle|u_\pm(k)\right\rangle
=\frac{i}{2}
$$

なので、

$$
\begin{aligned}
\gamma_\pm
&=-i\int_{-\pi}^{\pi}
\frac{\partial\phi}{\partial k}
\left\langle u_\pm(k)\middle|\frac{\partial}{\partial\phi}\middle|u_\pm(k)\right\rangle
\,\mathrm{d}k\\
&=\frac{1}{2}\int_{-\pi}^{\pi}
\frac{\partial\phi}{\partial k}\,\mathrm{d}k\\
&=\pi\nu\qquad (\mathrm{mod}\ 2\pi).
\end{aligned}
$$

したがって、このゲージと単位胞の取り方では \(\nu=1\) のとき \(\gamma_\pm=\pi\)、\(\nu=0\) のとき \(\gamma_\pm=0\) となる。Zak phase は一般には単位胞原点の選び方に依存するが、同じ規約で比較した相の差と、対称性による \(0/\pi\) の量子化が物理的に重要である。
"""

NUMERICAL_ZAK_SECTION = r"""### Zak phase の数値計算

数値固有ベクトルの位相は各 \(k\) で任意なので、Berry connection を有限差分して直接足し合わせる方法は不安定である。代わりに、閉じた離散波数列 \(k_0,\ldots,k_{N-1}\) に対して正規化したリンク変数

$$
U_{n,\pm}
=\frac{\langle u_\pm(k_n)\mid u_\pm(k_{n+1})\rangle}
{\left|\langle u_\pm(k_n)\mid u_\pm(k_{n+1})\rangle\right|},
\qquad k_N\equiv k_0
$$

を定義する。Wilson loop と Zak phase は

$$
W_\pm=\prod_{n=0}^{N-1}U_{n,\pm},
\qquad
\gamma_\pm=\operatorname{Arg}W_\pm
\quad (\mathrm{mod}\ 2\pi)
$$

で与えられる。リンクを絶対値で正規化することで、有限分割に由来する不要な振幅誤差を除き、位相だけを積算できる。最後の点と最初の点の重なりを必ず含める必要がある。
"""

GEN_DATA_CODE = r'''def gen_data(H, v, w, m=0.0, N=101):
    """Hermitian SSH model on N points in [-pi, pi)."""
    k_list = np.linspace(-np.pi, np.pi, N, endpoint=False)
    eigenvalues = []
    eigenvectors = []

    for k in k_list:
        Hamiltonian = H(k, v, w, m=m)
        E, U = np.linalg.eigh(Hamiltonian)
        eigenvalues.append(E)
        eigenvectors.append(U)

    eigenvectors = np.asarray(eigenvectors)
    return k_list, np.asarray(eigenvalues), eigenvectors, eigenvectors.copy()
'''

WILSON_CODE = r'''def zak_phase_from_eigenvectors(eigenvectors, atol=1e-12):
    """Return one Zak phase per band in the interval [-pi, pi]."""
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

zak_phase = zak_phase_from_eigenvectors(npr_psiR)
print("Zak phase / pi:", zak_phase / np.pi)
'''


def source_text(cell: dict) -> str:
    source = cell.get("source", [])
    return "".join(source) if isinstance(source, list) else str(source)


def set_source(cell: dict, text: str) -> None:
    cell["source"] = text.splitlines(keepends=True)


def replace_markdown_section(text: str, heading: str, replacement: str, keep_images: bool = False) -> str:
    pattern = re.compile(
        rf"(?ms)^{re.escape(heading)}\s*\n.*?(?=^#{{1,6}}\s|\Z)"
    )

    def repl(match: re.Match[str]) -> str:
        result = replacement.rstrip() + "\n"
        if keep_images:
            images = re.findall(r"!\[[^\]]*\]\([^)]+\)", match.group(0))
            if images:
                result += "\n" + "\n\n".join(dict.fromkeys(images)) + "\n"
        return result

    return pattern.sub(repl, text, count=1)


def normalize_markdown(text: str) -> str:
    # Commands unsupported by GitHub's math renderer.
    text = text.replace(r"\bm", r"\boldsymbol")
    text = re.sub(r"\\ket\{([^{}]+)\}", r"\\lvert \1 \\rangle{}", text)
    text = re.sub(r"\\bra\{([^{}]+)\}", r"\\langle \1 \\rvert{}", text)
    text = re.sub(r"\\rangle(?!\{\})", r"\\rangle{}", text)
    text = re.sub(r"\\rvert(?!\{\})", r"\\rvert{}", text)

    # Known factual and typographical corrections outside the rewritten sections.
    text = text.replace(
        r"$H_{1,2N}=H_{2N,1}=v$",
        r"$H_{1,2N}=H_{2N,1}=w$",
    )
    text = text.replace(
        "Trivialユニットセル内の2サイトの結合が強い状態、Non-trivialはその逆である。",
        "Trivial phase では単位胞内結合が強く、topological phase では単位胞間結合が強い（この判定は単位胞と終端の規約に依存する）。",
    )
    text = text.replace("Willson loop", "Wilson loop")
    text = text.replace("Hermite の場合", "Hermitian の場合")

    # GitHub Markdown does not render a bare equation environment.
    text = text.replace(r"\begin{equation}", "$$\n" + r"\begin{aligned}")
    text = text.replace(r"\end{equation}", r"\end{aligned}" + "\n$$")

    text = replace_markdown_section(text, "### 平面波基底による表現", PLANE_WAVE_SECTION)
    text = replace_markdown_section(text, "### 固有値と固有ベクトル", EIGENVECTOR_SECTION)
    text = replace_markdown_section(text, "### Chiral symmetry", CHIRAL_SECTION)
    text = replace_markdown_section(text, "### Winding number", WINDING_SECTION, keep_images=True)
    text = replace_markdown_section(text, "### Zak phase", ZAK_SECTION)
    text = replace_markdown_section(text, "### Zak phase の数値計算", NUMERICAL_ZAK_SECTION)
    return text


def clean_outputs(cell: dict) -> None:
    if cell.get("cell_type") != "code":
        return

    # The workflow does not execute the notebook. Remove stale outputs so the
    # generated README never shows results that no longer match the source.
    cell["outputs"] = []
    cell["execution_count"] = None


def normalize_code(text: str) -> str:
    if "def gen_data(" in text:
        return GEN_DATA_CODE
    if "wilson_loop = np.ones" in text:
        return WILSON_CODE

    text = text.replace(
        "plt.xticks(np.linspace(0,2*np.pi,5),[r'$-\\pi$',r'$-\\pi/2$',r'$0$',r'$\\pi/2$',r'$\\pi$'],fontsize=fs)",
        "plt.xticks(np.linspace(-np.pi, np.pi, 5), [r'$-\\pi$', r'$-\\pi/2$', r'$0$', r'$\\pi/2$', r'$\\pi$'], fontsize=fs)",
    )
    text = text.replace(
        "np.real(np.array(vecdat[vi][0::2,n].T))",
        "np.abs(np.asarray(vecdat[vi][0::2, n]).ravel())**2",
    )
    text = text.replace(
        "np.real(np.array(vecdat[vi][1::2,n].T))",
        "np.abs(np.asarray(vecdat[vi][1::2, n]).ravel())**2",
    )
    text = text.replace("plt.ylim(-1.2,1.2)", "plt.ylim(0, 1.0)")
    return text


def validate(notebook: dict) -> None:
    markdown = "\n".join(
        source_text(cell)
        for cell in notebook.get("cells", [])
        if cell.get("cell_type") == "markdown"
    )

    forbidden = [r"\\bm", r"\\ket", r"\\bra", r"\\begin\{equation\}"]
    for command in forbidden:
        if re.search(command, markdown):
            raise ValueError(f"Unsupported LaTeX command remains: {command}")

    if re.search(r"\\r(?:angle|vert)[A-Za-z]", markdown):
        raise ValueError("A closing Dirac-notation command is joined to following text.")
    if markdown.count("$$") % 2:
        raise ValueError("Unbalanced display-math delimiters ($$).")

    required = [
        r"\{\sigma_z,H(k)\}=0",
        r"\frac{1}{2}\int_{-\pi}^{\pi}",
        r"\gamma_\pm=\operatorname{Arg}W_\pm",
        r"H_{1,2N}=H_{2N,1}=w",
    ]
    for expression in required:
        if expression not in markdown:
            raise ValueError(f"Required corrected expression is missing: {expression}")

    forbidden_facts = [
        r"[\sigma_z, H(k)] = 0",
        "v<w$ なら $\\nu=0",
        r"H_{1,2N}=H_{2N,1}=v",
    ]
    for expression in forbidden_facts:
        if expression in markdown:
            raise ValueError(f"Known incorrect expression remains: {expression}")


def add_title_and_citation(notebook: dict) -> None:
    markdown_cells = [
        cell for cell in notebook.get("cells", [])
        if cell.get("cell_type") == "markdown"
    ]
    if not markdown_cells:
        raise ValueError("The notebook contains no Markdown cells.")

    first_cell = markdown_cells[0]
    first_text = source_text(first_cell)
    lines = first_text.splitlines()
    body_lines = []
    skip_reference_block = False

    for index, line in enumerate(lines):
        if index == 0 and line.startswith("# "):
            continue
        if line.strip() == "## Reference":
            skip_reference_block = True
            continue
        if skip_reference_block:
            if line.startswith("## ") and line.strip() != "## Reference":
                skip_reference_block = False
            else:
                continue
        if "arxiv.org/abs/1509.02295" in line:
            continue
        body_lines.append(line)

    body = "\n".join(body_lines).strip()
    set_source(first_cell, TITLE + "\n\n" + CITATION + body + "\n")


def main() -> None:
    if not NOTEBOOK_PATH.exists():
        raise FileNotFoundError(f"{NOTEBOOK_PATH} was not found.")

    with NOTEBOOK_PATH.open("r", encoding="utf-8") as f:
        notebook = json.load(f)

    for cell in notebook.get("cells", []):
        text = source_text(cell)
        if cell.get("cell_type") == "markdown":
            set_source(cell, normalize_markdown(text))
        elif cell.get("cell_type") == "code":
            set_source(cell, normalize_code(text))
            clean_outputs(cell)

    add_title_and_citation(notebook)
    validate(notebook)

    with NOTEBOOK_PATH.open("w", encoding="utf-8") as f:
        json.dump(notebook, f, ensure_ascii=False, indent=1)
        f.write("\n")


if __name__ == "__main__":
    main()
