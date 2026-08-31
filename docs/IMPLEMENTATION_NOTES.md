# SimplexLP — Implementation Notes

**Project #13 · Šehzada Sijarić (19964)**

This document maps the implementation to the proposal, derives the two
algorithms as implemented, records the design decisions and the validation
results.

---

## 1. Problem form and conversion

Everything is solved in standard form

```
min cᵀx   s.t.   A x = b,  x ≥ 0,     A ∈ ℝ^{m×n} sparse.
```

`LPProblem.h` converts the user's inequality model (`max/min`, rows with
`<= / >= / =`, implicit `x ≥ 0`):

* `maximize cᵀx` → `minimize (−c)ᵀx` (the sign is restored for display via
  `StandardLP::objSign`),
* `aᵀx ≤ b` → `aᵀx + s = b` (slack, coefficient **+1**),
* `aᵀx ≥ b` → `aᵀx − s = b` (surplus, coefficient −1).

`slackOfRow[i]` remembers which column is the +1 slack of row *i*: those
columns form an identity submatrix and provide a **crash basis** whenever
`b ≥ 0`, letting phase I be skipped — the common case for the classroom
`≤`-models in the 2D/3D tabs.

`SparseColMatrix` (LPCommon.h) stores A column-wise as `(row, value)` lists:
exactly the access pattern the revised simplex needs (columns `a_j` for
pricing and FTRAN) and exactly what `sparse::ISolver::addTriple` consumes.

## 2. Revised simplex (RevisedSimplex.h)

With basis index set `B` (one column per row position) and `B = A(:,B)`:

```
x_B = B⁻¹ b                  basic solution            (FTRAN)
y   = B⁻ᵀ c_B                simplex multipliers       (BTRAN)
d_j = c_j − yᵀ a_j           reduced costs             (sparse dots)
w   = B⁻¹ a_q                entering column direction (FTRAN)
θ   = min { x_B(i)/w_i : w_i > 0 }   ratio test → leaving row r
```

Only `B` is ever factorized — never the full tableau — which is the defining
property of the *revised* method and the reason it pairs naturally with a
sparse LU.

**natID mapping (BasisFactor.h).** Each (re)factorization creates two fresh
solvers via `sparse::createDblSolver(m, nnz, NonSymmetric, LU,
MarkowitzSinglePass)` — one loaded with `B`, one with `Bᵀ` — because
`ISolver` exposes `solveExt(rhs, x)` but no transpose solve and no
"reset values" call. `populateDiagonals(0.0)` seeds the diagonal pattern
first, mirroring the SDK's own MatrixTests; `addTriple` accumulates, so every
entry is inserted exactly once. This satisfies the proposal sentence "the
basis will be updated by solving the linear system B d = a_q ... using
natID's sparse LU solver" literally.

**Phase I with minimal artificials.** Rows are sign-normalized so `b ≥ 0`
(`rowSign` flips a row and its RHS; duals are mapped back at the end through
`y_orig = S y_w`). A row whose slack still has coefficient +1 starts with
that slack in the basis; only the remaining rows get an artificial with
phase-I cost 1. If no artificials are needed the phase-I loop never runs.
After phase I, `phase1Obj > tolFeas` ⇒ **Infeasible**; otherwise remaining
basic artificials are driven out with degenerate pivots: for the artificial
in row r, `z = B⁻ᵀ e_r` is one BTRAN, and `(B⁻¹ a_j)_r = zᵀ a_j` is a sparse
dot per candidate column. A row where no pivot exists is linearly dependent
(redundant); its artificial stays basic at 0 and is simply never priced in
phase II.

**Pivot rules.** Entering: Dantzig (most negative `d_j`); after `blandAfter`
consecutive degenerate pivots the code switches to **Bland's rule** (first
improving index), which provably terminates — verified on Beale's classic
cycling example in the test suite. Leaving: minimum ratio with ties broken by
the **largest |w_i|** for numerical stability.

**Recording for the visualization.** With `recordPath` on, every iteration
stores a snapshot `(iter, phase, entering, leaving, objective, x)`; the
canvases replay this sequence as the pivot animation, and the IPM stores its
iterates the same way, so both trajectories are drawn from identical data
structures.

**Complexity note.** One LU per pivot (plus one for `Bᵀ`) is the clear,
robust choice and showcases the natID solver; classic production codes amortize
with eta-files / Forrest–Tomlin updates and periodic refactorization — listed
under future work.

## 3. Interior-point method (InteriorPoint.h)

Primal–dual pair:

```
(P) min cᵀx, Ax = b, x ≥ 0        (D) max bᵀy, Aᵀy + s = c, s ≥ 0
```

Mehrotra predictor–corrector on the **normal equations**. With
`D² = diag(x/s)` and residuals `r_p = b − Ax`, `r_d = c − Aᵀy − s`:

```
(A D² Aᵀ) Δy = r_p + A S⁻¹ (X r_d − r_xs)
Δs = r_d − Aᵀ Δy
Δx = S⁻¹ (−X Δs − r_xs)
```

The affine pass uses `r_xs = XSe`; the corrector reuses the same
factorization with `r_xs = XSe − σμe + ΔX_aff ΔS_aff e`,
`σ = (μ_aff/μ)³` — one factorization, two solves per iteration.
Steps take `0.99×` the fraction-to-boundary; the start point is Mehrotra's
(two solves with `AAᵀ`, then a positivity shift).

**natID mapping.** `M = A D² Aᵀ` is symmetric positive definite, so the
default solver is created with `Symmetry::SymmetricPosDef` and **only the
lower triangle inserted** — the convention confirmed in the SDK's
`MatrixTests`. The factorization engine is `SolverType::LU`: in this SDK
build that is the supported path for SPD systems (the SDK's own Cholesky
test, `TestCholesky.h`, likewise creates SPD solvers with `SolverType::LU`),
whereas `SolverType::LLT` aborts inside `createDblSolver`. The sparsity
pattern of `M` is computed once (symbolically, from column outer products
`Σ_j d_j a_j a_jᵀ`) and only the values are refilled each iteration; a tiny
diagonal regularization (`1e-10`-scaled) keeps the factorization stable, and
if the SPD factorization ever fails the code transparently falls back to
plain LU with mirrored entries.

The method assumes feasible-and-bounded problems — guaranteed by the
generator's construction; a homogeneous self-dual embedding would extend it
to infeasibility detection (future work).

## 4. Synthetic generator (LPGenerator.h)

Instances are built with **both certificates**, so an optimum always exists:

* sparse `A`: a round-robin anchor guarantees every row is touched, the rest
  of each column is filled randomly (`nnzPerCol` target, values U(−1,1));
* primal: pick `x* > 0`, set `b = A x*`  ⇒ primal feasible;
* dual: pick `y*` free, `s* > 0`, set `c = Aᵀy* + s*`  ⇒ dual feasible.

By strong duality the instance is feasible and bounded. The optimal value is
not `cᵀx*` (the simplex typically finds a better vertex) — correctness is
checked through KKT residuals and the simplex-vs-IPM objective gap instead.

## 5. Geometry for the visualization (Polytope.h)

* **2D:** the feasible polygon is obtained by Sutherland–Hodgman clipping of
  a large box against every half-plane (`=` contributes both directions).
* **3D:** vertex enumeration over all plane triples (user planes + axis
  planes + a safety box that caps unbounded regions), feasibility filtering,
  deduplication; each plane's face is its incident vertices sorted by angle
  in the plane basis. The canvas renders faces back-to-front (painter's
  algorithm) with the model rotated by yaw/pitch.
* The 2×2 / 3×3 intersection systems are solved with **`dense::Matrix`**
  (`getManipulator`, `solve`) — the proposal's dense-matrix usage, in the
  place where a dense solve is genuinely the right tool.

## 6. MPS and the CLP comparison

`MPS.h` writes `NAME / ROWS (N,E,L,G) / COLUMNS / RHS / ENDATA` in the
fixed-field layout CLP accepts and reads the same subset back (slack/surplus
conversion on input; RANGES/BOUNDS are rejected explicitly). The benchmark
(GUI tab or `lpBench`) exports every generated instance, and
`scripts/run_clp.sh` runs CLP twice per file (`-primalsimplex`, `-barrier`),
collecting time and objective into `clp_results.csv` — so the natID
implementation and CLP are measured on **identical inputs**, as the proposal
requires.

## 7. GUI architecture

```
MainWindow ── MenuBar / ToolBar
    └── MainView (StandardTabView)
          ├── Tab2D  = SplitterLayout( ViewVizControls │ View2DCanvas )
          ├── Tab3D  = SplitterLayout( ViewVizControls │ View3DCanvas )
          └── TabBench = SplitterLayout( parameters │ 2 × ViewChartCanvas )
```

* `VizModel` owns the parsed `UserModel`, the `StandardLP`, both `LPResult`s
  (paths recorded), the polygon/polytope geometry and the animation cursor;
  the canvases only read it.
* `LinearExprParser` accepts `3x + 2y <= 18`-style lines (variables
  `x,y,z` or `x1..x3`, both-side expressions, comments).
* Animation uses `gui::Timer` on the tab view (`onTimer` advances the
  cursor); the Run menu and the toolbar route through
  `MainWindow::onActionItem` → `MainView` → active tab.
* The benchmark sweep runs on a `std::thread`; progress is marshalled with
  `gui::thread::asyncExecInMainThread` (the SDK's animation-example pattern),
  the UI only ever touches widgets on the main thread.
* Settings dialog = the SDK's language-switcher pattern (EN/BA translations
  in `res/tr/`, restart prompt on change).

**Chart widget decision.** The benchmark charts use a small custom
`ViewChartCanvas` (axes with "nice" ticks, multi-series polylines, legend)
instead of natPlot: the uploaded SDK build was missing `gui/plot/Plot.h` and
no shipped example exercises `gui::plot::View`, so its API could not be
verified end-to-end. The custom widget uses only Canvas/Shape/DrawableString
calls that the SDK examples themselves use. Swapping in
`gui::plot::View(fontAxis, fontLegend, colorScheme)` +
`addFunction(x, y, len, color, ...)` later is a localized change in
`TabBench`.

## 8. Validation

`tests/standalone/` compiles the **unmodified** core headers against shim
implementations of `sparse::ISolver` (dense LU with the same accumulate /
populateDiagonals / one-triangle-for-symmetric semantics), `dense::Matrix`
and `cnt::SafeFullVector`, then runs **86 checks**:

* Dantzig's production example: optimum 36 at (2,6), crash basis ⇒ 0 phase-I
  iterations, objective monotone non-decreasing along the path;
* a `≥`-model exercising phase I; an equality model; explicit infeasible and
  unbounded detection; a degenerate vertex; **Beale's cycling example**
  (terminates at −0.05 thanks to the Bland fallback); duplicated equality
  rows (redundancy handling);
* MPS write→read round trip reproduces the instance exactly;
* 2D polygon and 3D cube/tetrahedron geometry (exact vertex/face counts);
* 8 random generated instances up to (m=60, n=150): simplex and IPM both
  `Optimal`, relative objective gaps 1e-8…1e-10, KKT residuals pass.

Headline comparison on those instances (shim LU, so absolute times are not
meaningful — iteration counts are): simplex grows from 10 to 182 iterations
while the IPM stays at 7–11, the textbook contrast the benchmark tab then
shows at scale with the real natID factorizations.

Additionally, the complete GUI application and the console tool were
syntax-checked (`g++ -fsyntax-only -std=c++20`) against the real SDK
headers, so every natID call site (signatures, enums, access levels) is
verified even before linking on the target machine.

## 9. Future work

* Eta-file / Forrest–Tomlin basis updates with periodic refactorization;
* steepest-edge or Devex pricing;
* Harris two-pass ratio test;
* homogeneous self-dual IPM (infeasibility certificates);
* presolve (empty/singleton rows & columns, bound tightening);
* full MPS BOUNDS/RANGES support for Netlib problems.
