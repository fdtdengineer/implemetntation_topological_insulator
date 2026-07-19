# 1次元トポロジカル絶縁体：SSH模型の実装

## Reference

J. K. Asbóth, L. Oroszlány, and A. Pályi, *A Short Course on Topological Insulators: Band-Structure Topology and Edge States in One and Two Dimensions*, Lecture Notes in Physics **919** (Springer, 2016). [doi:10.1007/978-3-319-25607-8](https://doi.org/10.1007/978-3-319-25607-8); [arXiv:1509.02295](https://arxiv.org/abs/1509.02295).

## 概要

Su–Schrieffer–Heeger（SSH）模型は、単位胞ごとに A・B の2サイトを持つ1次元の tight-binding model である。単位胞内結合を $v$、単位胞間結合を $w$ とする。

![SSH chain](README_files/image.png)

以下では、標準的な単位胞、格子定数 $a=1$、実数結合 $v,w>0$ を基準とする。

## 実空間ハミルトニアン

$N$ 個の単位胞からなる開放鎖のハミルトニアンは

$$
\begin{aligned}
H={}&v\sum_{l=1}^{N}
\left(
|l,B\rangle\langle l,A|+|l,A\rangle\langle l,B|
\right)\\
&+w\sum_{l=1}^{N-1}
\left(
|l+1,A\rangle\langle l,B|+|l,B\rangle\langle l+1,A|
\right).
\end{aligned}
$$

基底を

$$
\left(|1,A\rangle,|1,B\rangle,|2,A\rangle,|2,B\rangle,\ldots\right)
$$

と並べると、$N=4$ の行列は

$$
H=
\begin{bmatrix}
0&v&0&0&0&0&0&0\\
v&0&w&0&0&0&0&0\\
0&w&0&v&0&0&0&0\\
0&0&v&0&w&0&0&0\\
0&0&0&w&0&v&0&0\\
0&0&0&0&v&0&w&0\\
0&0&0&0&0&w&0&v\\
0&0&0&0&0&0&v&0
\end{bmatrix}.
$$

周期境界条件では、最後の B サイトと最初の A サイトを結ぶため

$$
H_{1,2N}=H_{2N,1}=w
$$

を追加する。

## Bloch Hamiltonian

周期境界条件の下で

$$
|k,\alpha\rangle
=\frac{1}{\sqrt{N}}\sum_{l=1}^{N}e^{ikl}|l,\alpha\rangle,
\qquad \alpha\in\{A,B\}
$$

と定義する。基底を $\left(|k,A\rangle,|k,B\rangle\right)$ とすると

$$
\begin{aligned}
H(k)
&=
\begin{bmatrix}
0&v+we^{-ik}\\
v+we^{ik}&0
\end{bmatrix}\\
&=
\begin{bmatrix}
0&h(k)^*\\
h(k)&0
\end{bmatrix},
\qquad h(k)=v+we^{ik}.
\end{aligned}
$$

したがって

$$
H(k)=d_x(k)\sigma_x+d_y(k)\sigma_y,
\qquad
d_x(k)=v+w\cos k,\quad d_y(k)=w\sin k.
$$

バンドギャップは $h(k)=0$ で閉じる。$v,w>0$ なら転移点は

$$
v=w,\qquad k=\pi
$$

である。

## 固有値と固有ベクトル

$h(k)=|h(k)|e^{i\phi(k)}$ と書くと

$$
E_\pm(k)
=\pm|h(k)|
=\pm\sqrt{v^2+w^2+2vw\cos k},
$$

$$
|u_\pm(k)\rangle
=\frac{1}{\sqrt{2}}
\begin{bmatrix}
1\\
\pm e^{i\phi(k)}
\end{bmatrix}.
$$

位相は

$$
\phi(k)=\mathrm{Arg}\!\left(v+we^{ik}\right)
=\mathrm{atan2}\!\left(w\sin k,\,v+w\cos k\right)
$$

で求める。単純な $\arctan(y/x)$ では象限と分岐を正しく追跡できない。ギャップ閉鎖点 $h(k)=0$ では、この位相表示は定義できない。

## カイラル対称性

$H(k)$ に $\sigma_z$ 成分がないため

$$
\{\sigma_z,H(k)\}=0,
$$

すなわち

$$
\sigma_zH(k)\sigma_z^{-1}=-H(k)
$$

が成り立つ。これは**反交換関係**であり、$[\sigma_z,H(k)]=0$ ではない。

また、$\sigma_z=\mathrm{diag}(1,-1)$ は A/B サブ格子へ相対的な符号を与える演算子である。A/B を交換する演算子は $\sigma_x$ である。

## Winding number

ギャップが開いている $h(k)\neq0$ の場合、winding number を

$$
\nu
=\frac{1}{2\pi i}\int_{-\pi}^{\pi}
\frac{d}{dk}\log h(k)\,dk
=\frac{1}{2\pi}\int_{-\pi}^{\pi}
\frac{\partial\phi(k)}{\partial k}\,dk
$$

と定義する。これは $h(k)$ の軌跡が複素平面の原点を向き付きで何周するかを表す。

この Fourier 規約と $v,w>0$ の下では

$$
\begin{aligned}
\nu&=0 && (v>w,\ \mathrm{trivial}),\\
\nu&=1 && (v<w,\ \mathrm{topological}).
\end{aligned}
$$

$v=w$ ではギャップが閉じるため、winding number は定義できない。実数結合の符号も許す場合、この規約では $|w|>|v|$ なら $\nu=1$、$|w|<|v|$ なら $\nu=0$ である。Fourier変換の向きを逆にすると $\nu$ の符号は反転するが、相の分類は変わらない。

![Winding trajectory](README_files/image-4.png)

## Zak phase

Berry connection を

$$
A_\pm(k)
=-i\left\langle u_\pm(k)\right|
\frac{\partial}{\partial k}
\left|u_\pm(k)\right\rangle
$$

と定義し、Zak phase を

$$
\gamma_\pm
=\int_{-\pi}^{\pi}A_\pm(k)\,dk
\qquad (\mathrm{mod}\ 2\pi)
$$

とする。上の固有ベクトルに対して

$$
\left\langle u_\pm(k)\right|
\frac{\partial}{\partial\phi}
\left|u_\pm(k)\right\rangle
=\frac{i}{2}
$$

なので

$$
\begin{aligned}
\gamma_\pm
&=-i\int_{-\pi}^{\pi}
\frac{\partial\phi}{\partial k}
\left\langle u_\pm(k)\right|
\frac{\partial}{\partial\phi}
\left|u_\pm(k)\right\rangle dk\\
&=\frac{1}{2}\int_{-\pi}^{\pi}
\frac{\partial\phi}{\partial k}\,dk\\
&=\pi\nu
\qquad (\mathrm{mod}\ 2\pi).
\end{aligned}
$$

したがって、このゲージと単位胞の規約では

$$
\nu=0\Rightarrow\gamma_\pm=0,
\qquad
\nu=1\Rightarrow\gamma_\pm=\pi
\quad (\mathrm{mod}\ 2\pi).
$$

Zak phase 自体は単位胞原点の選び方に依存する。物理的に重要なのは、同じ規約で比較した位相差と、対称性による $0/\pi$ の量子化である。

## バルク–エッジ対応

異なる winding number を持つ2領域の境界には、カイラル対称性が保たれている限り

$$
|\Delta\nu|
$$

個の零エネルギー境界モードが現れる。有限の topological SSH 鎖では左右の端に1個ずつ端状態が現れ、有限サイズでは両者の混成によりエネルギーが指数関数的にわずかに分裂する。

端状態の有無は、bulk の単位胞規約と実際の鎖の終端を整合させて判断する必要がある。

## Python実装

### Bloch Hamiltonian

```python
import numpy as np
import matplotlib.pyplot as plt


def H_SSH(k: float, v: float, w: float, m: float = 0.0) -> np.ndarray:
    """Bloch Hamiltonian in the (A, B) sublattice basis."""
    return np.array(
        [
            [m, v + w * np.exp(-1j * k)],
            [v + w * np.exp(1j * k), -m],
        ],
        dtype=complex,
    )
```

$m\sigma_z$ は staggered onsite term であり、$m\neq0$ ではカイラル対称性を破る。この場合、Zak phase は一般に $0$ または $\pi$ に量子化されず、零エネルギー端状態もカイラル対称性によって保護されない。

### バンド構造

```python
def band_data(v: float, w: float, m: float = 0.0, n_k: int = 401):
    k_grid = np.linspace(-np.pi, np.pi, n_k, endpoint=False)
    eigenvalues = []
    eigenvectors = []

    for k in k_grid:
        values, vectors = np.linalg.eigh(H_SSH(k, v, w, m))
        eigenvalues.append(values)
        eigenvectors.append(vectors)

    return k_grid, np.asarray(eigenvalues), np.asarray(eigenvectors)


k_grid, eigenvalues, eigenvectors = band_data(v=0.8, w=1.0)

for band in range(2):
    plt.plot(k_grid, eigenvalues[:, band])
plt.xlabel(r"$k$")
plt.ylabel(r"$E$")
plt.xticks(
    [-np.pi, -np.pi / 2, 0, np.pi / 2, np.pi],
    [r"$-\pi$", r"$-\pi/2$", r"$0$", r"$\pi/2$", r"$\pi$"],
)
plt.show()
```

### Gauge-invariant な Zak phase

数値固有ベクトルの位相は各 $k$ で任意なので、Berry connection を単純に有限差分して足し合わせる方法は不安定である。代わりに

$$
U_{n,\pm}
=\frac{\langle u_\pm(k_n)|u_\pm(k_{n+1})\rangle}
{|\langle u_\pm(k_n)|u_\pm(k_{n+1})\rangle|},
\qquad k_N\equiv k_0
$$

を用いる。Wilson loop と Zak phase は

$$
W_\pm=\prod_{n=0}^{N-1}U_{n,\pm},
\qquad
\gamma_\pm=\mathrm{Arg}\,W_\pm
$$

である。

```python
def zak_phase_from_eigenvectors(
    eigenvectors: np.ndarray,
    atol: float = 1e-12,
) -> np.ndarray:
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


for v, w in [(1.2, 1.0), (0.8, 1.0)]:
    _, _, vectors = band_data(v=v, w=w)
    print(v, w, zak_phase_from_eigenvectors(vectors) / np.pi)
```

期待される結果は、trivial相で $0$、topological相で $\pm1$、すなわち $\pm\pi\equiv\pi\pmod{2\pi}$ である。

### 開放有限鎖

```python
def H_SSH_open(n_cells: int, v: float, w: float) -> np.ndarray:
    """Open SSH chain in the basis (A1, B1, A2, B2, ...)."""
    H = np.zeros((2 * n_cells, 2 * n_cells), dtype=float)

    for cell in range(n_cells):
        a, b = 2 * cell, 2 * cell + 1
        H[a, b] = H[b, a] = v

        if cell < n_cells - 1:
            next_a = 2 * (cell + 1)
            H[b, next_a] = H[next_a, b] = w

    return H


n_cells = 24
for v, w in [(1.5, 1.0), (0.5, 1.0)]:
    values, vectors = np.linalg.eigh(H_SSH_open(n_cells, v, w))
    near_zero = np.count_nonzero(np.abs(values) < 1e-6)
    print(v, w, near_zero)

    edge_index = np.argmin(np.abs(values))
    probability = np.abs(vectors[:, edge_index]) ** 2
    plt.plot(probability, marker="o")
    plt.xlabel("site index")
    plt.ylabel(r"$|\psi|^2$")
    plt.show()
```

十分長い鎖では、topological相 $v<w$ に2個の near-zero mode が現れ、trivial相 $v>w$ には現れない。

## 自動検証

`scripts/verify_ssh.py` では、次の3段階を自動確認する。

1. Hermiticity、解析固有値、カイラル反交換関係。
2. trivial/topological両相のWilson-loop Zak phase。
3. 開放有限鎖における端状態の有無。

README生成時にも、GitHubで未対応の数式マクロ、インライン数式区切り、既知の誤式が残っていないかを検査する。

## 参考文献

- J. K. Asbóth, L. Oroszlány, and A. Pályi, *A Short Course on Topological Insulators*, arXiv:1509.02295.
- T. Fukui, Y. Hatsugai, and H. Suzuki, “Chern Numbers in Discretized Brillouin Zone,” *J. Phys. Soc. Jpn.* **74**, 1674 (2005).
