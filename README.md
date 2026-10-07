# Static response, structure factor and pair correlation of the 3D electron gas from DMC

Data accompanying the joint submission:

- M. H. Güneş, Y. Yang, V. Gorelov, M. A. Morales, M. Gatti, L. Reining, S. Zhang,
  "Static Response of the Strongly Coupled Electron Gas from Diffusion Monte Carlo", submitted to Phys. Rev. B.
- M. H. Güneş, Y. Yang, V. Gorelov, M. A. Morales, M. Gatti, L. Reining, S. Zhang,
  "Signature of Short-Range Order in Static Response of the Three-Dimensional Electron Gas",
  submitted to Phys. Rev. Lett., [arXiv:2609.10006](https://arxiv.org/abs/2609.10006).

Citation metadata: `CITATION.cff`.

All files are plain text, whitespace-separated, one `#` header line naming the columns.
Hartree atomic units throughout.

## Conventions

- rs: Wigner–Seitz radius (bohr); unpolarized 3D homogeneous electron gas, simple cubic simulation cell.
- n0 = 3 / (4π rs³) — electron density.
- kF = (3π² n0)^(1/3) = (9π/4)^(1/3) / rs — Fermi wavevector.
- Allowed wavevectors in a cell of Ne electrons: q = (2π/L) n, n ∈ ℤ³, L = (4π rs³ Ne / 3)^(1/3).

## Files

| files | content | columns |
|---|---|---|
| `G/G_rs{rs}.txt` | static local field factor G(q) | `q/kF  G  dG  Ne` |
| `chi/chi_rs{rs}.txt` | static density response χ(q) | `q/kF  -chi/n0  dchi/n0  Ne` |
| `sk/sk_rs{rs}_Ne{Ne}.txt` | static structure factor S(q) | `q/kF  S  dS` |
| `gofr/gofr_rs{rs}_Ne{Ne}.txt` | pair correlation function g(r) | `r*kF  g  dg` |

Available rs values:

| quantity | rs | Ne |
|---|---|---|
| G, χ | 2, 5, 10, 70, 90 | 162 |
| G, χ | 20, 50, 100 | 54 |
| G, χ | 30 | mixed: 16 points with Ne = 54, 6 points with Ne = 162 (q/kF ≈ 3.0–3.9) |
| S | 20, 30, 50, 80, 90 | 162 |
| S | 80, 90, 100, 110 | 250 |
| g | 20, 30, 50, 80, 90 | 162 |
| g | 5 | 54 |

For G and χ, Ne is given per row (`Ne` column). Since L(Ne=162)/L(Ne=54) = 3^(1/3), no q value is shared by the two sizes.

## G(q) and χ(q)

These are exactly the DMC points shown in the paper's figures.

- χ(q) is the DMC static density response, extracted from DMC energies E(v_q) at several perturbation
  amplitudes v_q and finite-size corrected via
  χ⁻¹(q,∞) = χ⁻¹(q,N) + [χ0⁻¹(q,∞) − χ0⁻¹(q,N)] (see paper).
- G(q) is defined by the Dyson equation
  χ(q) = χ0(q) / [1 − Vc(q)(1 − G(q)) χ0(q)], Vc(q) = 4π/q²,
  i.e. G(q) = 1 + [1/χ(q) − 1/χ0(q)] / Vc(q),
  with χ0 the Lindhard function (both spins),
  χ0(q) = −(kF/π²) [1/2 + (1 − x²)/(4x) ln|(1 + x)/(1 − x)|], x = q/(2kF).
- `chi_rs*.txt` gives −χ(q)/n0 (units Ha⁻¹; positive), so χ(q) = −n0 × (column 2).
- Errors (1σ): dχ from a parametric bootstrap — each E_DMC(v_q) is resampled within its error bar,
  E(v_q) is refit and the finite-size correction applied; dχ is the standard deviation over samples.
  dG = dχ / (Vc χ²).
- Trial wavefunction: Slater–Jastrow (SJ) for all points, except the five largest-q points at rs = 2
  (q/kF ≈ 1.90–2.58), which use Slater–Jastrow + backflow (SJ+BF).

## S(q) and g(r)

- Extrapolated estimators: S = 2 S_DMC − S_VMC, g = 2 g_DMC − g_VMC, where S_DMC, g_DMC are DMC mixed estimators
  and S_VMC, g_VMC the corresponding VMC estimators; this linear extrapolation removes the leading-order
  trial-wavefunction bias of the mixed estimator.
- Errors (1σ): dS = √(4 dS_DMC² + dS_VMC²), dg = √(4 dg_DMC² + dg_VMC²) (DMC and VMC errors treated as independent).
- S(q) is given at the allowed wavevectors of the Ne-electron cell.
- g(r) is histogrammed in uniform bins up to r ≈ L/2; `r*kF` is the bin center.
  At large rs, the smallest-r bins can be empty (g = 0, dg = 0), and a few small-r bins at rs = 80, 90
  are slightly negative (|g| ~ 10⁻⁵) due to the extrapolation; they are kept as is.

## Reading

```python
import numpy as np
x, G, dG, Ne = np.loadtxt("G/G_rs20.txt", unpack=True)
x, S, dS = np.loadtxt("sk/sk_rs80_Ne250.txt", unpack=True)
```

## License

This data set is licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) (see `LICENSE`).
