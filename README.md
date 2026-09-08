# SimplexLP

An implementation of the revised simplex method and a primal dual interior
point method for linear programming, both built on the natID SDK's sparse LU
solver, with a GUI that draws what the two methods actually do.

Course project #13, Šehzada Sijarić (19964).

Feed it a small model such as

```
max 3x + 5y
x <= 4
2y <= 12
3x + 2y <= 18
```

and the 2D tab draws the feasible polygon and walks the simplex from vertex to
vertex, with the interior point trajectory optionally overlaid cutting straight
through the middle. That contrast, boundary against interior, is the reason the
visualisation exists.

## What it does

Three tabs.

The 2D and 3D tabs take a model typed as ordinary text, one constraint per
line, and solve it with both methods while recording every iterate. The
feasible region is computed geometrically: Sutherland Hodgman clipping in two
dimensions, vertex enumeration over plane triples with painter's algorithm
rendering in three. Playback can be stepped or animated.

The benchmark tab sweeps synthetic instances of growing size, solves each with
both methods, and charts time and iteration counts. The same sweep is available
as a console tool, `lpBench`, which also exports every instance as MPS so
`scripts/run_clp.sh` can run COIN-OR CLP over identical inputs for an external
comparison.

## Implementation

The revised simplex keeps only the basis factorized, never the full tableau,
which is what makes a sparse LU the natural fit. Each refactorization builds
two natID solvers, one for the basis and one for its transpose, since the SDK
exposes no transpose solve. Phase I introduces artificials only on rows that
have no usable slack, and is skipped entirely when the slack columns already
give a crash basis. Pricing uses Dantzig's rule, falling back to Bland's rule
after a run of degenerate pivots so that termination is guaranteed, which the
test suite confirms on Beale's cycling example.

The interior point method is Mehrotra predictor corrector on the normal
equations. The matrix `A D² Aᵀ` is symmetric positive definite, so it is
handed to natID as a lower triangle only, its sparsity pattern computed once
and refilled each iteration. One factorization serves both the affine and the
corrector solve.

Full derivations, the natID API mapping and the design decisions behind both
are in [docs/IMPLEMENTATION_NOTES.md](docs/IMPLEMENTATION_NOTES.md).

## Validation

`tests/standalone/` compiles the unmodified solver headers against shim
implementations of the SDK interfaces and runs 86 checks: Dantzig's production
example, phase I and equality models, infeasible and unbounded detection,
degeneracy, Beale cycling, redundant rows, MPS round tripping, exact vertex and
face counts for 2D and 3D geometry, and eight random instances up to m=60,
n=150 where simplex and interior point agree to within 1e-8.

On those instances the simplex iteration count climbs from 10 to 182 as the
problems grow while the interior point method stays between 7 and 11, which is
the textbook contrast the benchmark tab then reproduces at scale.

```bash
cd tests/standalone && ./build_and_run.sh
```

## Installing

Prebuilt installers for Windows, macOS and Linux are on the
[releases page](https://github.com/ssijaric1/SimplexLP/releases).

| Platform | File | Install |
|---|---|---|
| Windows | `SimplexLP-win.zip` | unzip, run the `.exe` (keep the `.msi` beside it) |
| macOS | `SimplexLP-macOS-Silicon.zip` or `SimplexLP-macOS-Intel.zip` | unzip, drag to Applications, first launch needs right click then Open |
| Linux | `SimplexLP-linux.zip` | unzip, then `sudo apt install ./simplexLP.deb` |

Use `apt install` rather than `dpkg -i` on Linux, since dpkg will not pull in
the GTK dependencies. The macOS bundle is unsigned, which is why the first
launch needs the right click.

## Building from source

Needs the natID SDK at `~/natID.SDK` and GTK 4.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

This produces both `simplexLP` (the GUI) and `lpBench` (the console benchmark)
in `$RAMDisk/Out/simplexLPSol/Release`.

## Layout

```
src/core/     RevisedSimplex, InteriorPoint, BasisFactor, LPProblem, MPS,
              LPGenerator, Polytope, Verify
src/viz/      2D and 3D tabs, canvases, model parser
src/bench/    benchmark tab, sweep runner, chart canvas
cli/          lpBench console tool
scripts/      COIN-OR CLP comparison
tests/        standalone test suite with SDK shims
installer/    SetupCollector config and the patched GTK4 manifest
```

The interface is translated into English and Bosnian.

## Licence

MIT, see [LICENSE](LICENSE).
