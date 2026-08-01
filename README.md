# SeBa — `cferrus` fork

SeBa is a software package for simulating the evolution of single and binary
stars from the zero-age main sequence through to remnant phases. It is valid for
masses in the range 0.01–100 M☉ with variable metallicity, and includes
prescriptions for stellar-wind mass loss, supernovae, and binary interactions
(mass transfer, common envelope, tidal circularisation, gravitational-wave
inspiral, and magnetic braking).

This is a **modified fork** of [`amusecode/SeBa`](https://github.com/amusecode/SeBa),
tuned for large population-synthesis runs. It differs from upstream in several
important ways — see [What is different in this fork](#what-is-different-in-this-fork).

> For deeper implementation notes (build internals, thread-safety scope,
> optimisation details, accuracy benchmarks, cluster workflow) see
> [`CLAUDE.md`](CLAUDE.md).

## Contents

1. [What is different in this fork](#what-is-different-in-this-fork)
2. [Build](#build)
3. [Running SeBa](#running-seba)
4. [Output format](#output-format)
5. [Filtering the output](#filtering-the-output)
6. [Performance](#performance)
7. [Known caveats](#known-caveats)
8. [References](#references)

---

## What is different in this fork

| Area | Upstream `amusecode/SeBa` | This fork |
|------|---------------------------|-----------|
| **Output columns** | 18 columns per line | **9 columns** per line (see below) |
| **Lines per binary** | one per evolutionary event | **exactly 2**: initial state (T=0) and final state (T=end) |
| **Terminal output** | verbose stderr/stdout on every phase transition | **silent by default**; opt in with `-V` |
| **`binev.data`** | written | never created |
| **Parallelism** | serial | **OpenMP** parallel binary loop (Linux/GCC) |
| **File I/O** | opens/closes output file per binary | opens output once, buffers per-thread |
| **Hot-path math** | `pow()` everywhere | `cbrt()` / cached constants / multiplies |

The net effect is a smaller, faster, quieter tool that emits a compact 9-column
snapshot suitable for downstream population analysis, rather than a full
event-by-event history.

---

## Build

SeBa requires only a C++ compiler and `make`. No external dependencies.

```bash
# From the repo root:
make clean && make        # builds libsstar.a, libnode.a, libstd.a, librdc.a
cd dstar && make          # builds libdstar.a and the SeBa executable
```

The executable lands at `dstar/SeBa`.

- **Linux / GCC**: the `dstar/Makefile` detects `uname -s == Linux` and adds
  `-fopenmp` automatically, enabling multi-core parallelism.
- **macOS / Apple Clang**: OpenMP is not available without `libomp`; the build
  falls back to correct single-threaded execution automatically.
- **macOS / Homebrew GCC**: force OpenMP with `CXX=g++-14 make`.

### Relinking after library changes

`dstar/SeBa` only relinks when `SeBa.C` is newer than the binary. After editing
any library file (e.g. `double_star.C`, `single_star.C`), force a relink:

```bash
rm -f dstar/SeBa && cd dstar && make SeBa
```

### Smoke test

```bash
cd dstar
./SeBa -M 2 -m 1 -e 0.2 -a 200 -T 13500 -z 0.001
wc -l SeBa.data   # must be 2
```

---

## Running SeBa

Run `./SeBa -h` for the full CLI reference (all distribution options and
defaults). The three main modes:

### Single system

Primary mass M=2 M☉, secondary m=1 M☉, eccentricity e=0.2, separation a=200 R☉,
time T=13500 Myr, metallicity z=0.001:

```bash
./SeBa -M 2 -m 1 -e 0.2 -a 200 -T 13500 -z 0.001
```

### Input file

One binary per line as `a  e  M  m  z` (separation, eccentricity, primary mass,
secondary mass, metallicity):

```bash
./SeBa -I /absolute/path/to/input.txt -T 12550
```

```
200 0.2 2   1   0.001
500 0.5 2.5 1.5 0.02
```

> **Always use an absolute path for `-I`.** A path that cannot be opened causes a
> silent infinite loop in `read_binary_params` (a CPU core pegs indefinitely with
> no error message).

### Random population (Monte Carlo)

```bash
./SeBa -R -n 250000 -T 13500 -f 4 -m 0.95 -M 10
```

Initial parameters are drawn from selectable distributions (mass function,
semi-major axis, eccentricity, mass ratio). See `./SeBa -h`.

### Key flags

| Flag | Meaning |
|------|---------|
| `-h` | Print full help and exit |
| `-I file` | Read initial conditions from file (use absolute path) |
| `-R` | Generate random initial conditions |
| `-n N` | Number of binaries (random mode) |
| `-N offset` | Starting binary identity number |
| `-T Myr` | End time in Myr |
| `-s seed` | Fix random seed (required for reproducibility) |
| `-z Z` | Metallicity (0.0001–0.03) |
| `-O file` | Output filename (default: `SeBa.data`) |
| `-D` | Stop at merger or disruption |
| `-S` | Stop at first remnant formation |
| `-V` | Verbose: print diagnostics to stderr/stdout |

### Verbose mode

By default SeBa produces **no terminal output**. Add `-V` to enable diagnostics:

```bash
./SeBa -V -M 2 -m 1 -e 0.2 -a 200 -T 13500 -z 0.001
```

With `-V`, stdout carries the run log (seed, CE/kick settings, distribution
functions) and stderr carries per-component phase transitions, merger/spiral-in
events, kick diagnostics, and physics warnings. Omit `-V` for all batch runs.

### Reproducibility under OpenMP

The RNG (used for NS/BH kick angles) is not thread-safe. For a reproducible
seeded run, pin to a single thread:

```bash
OMP_NUM_THREADS=1 ./SeBa -I /abs/input.txt -T 12550 -s 42
```

---

## Output format

`SeBa.data` is written with **9 columns**. Each binary produces exactly **2
lines**: the initial state (T=0) and the final state (T=end). No intermediate
timesteps are written.

| Col | Field | Units |
|-----|-------|-------|
| 1 | binary identity | integer |
| 2 | binary type | integer (enum) |
| 3 | time | Myr |
| 4 | semi-major axis | R☉ |
| 5 | eccentricity | — |
| 6 | primary stellar type | integer |
| 7 | primary mass | M☉ |
| 8 | secondary stellar type | integer |
| 9 | secondary mass | M☉ |

### Binary type integers

| Value | Meaning |
|-------|---------|
| 0 | Unset |
| 2 | Detached |
| 3 | Semi-detached (stable mass transfer) |
| 7 | Disrupted |
| 9 | Merged |

### Stellar type integers

| Value | Type |
|-------|------|
| 3 | Main sequence |
| 5 | Hertzsprung gap |
| 6 | Sub-giant (first-ascent red giant) |
| 7 | Core-helium-burning (horizontal branch) |
| 8 | AGB / super-giant |
| 10 | Helium star |
| 11 | Helium giant |
| 12–14 | White dwarf (CO / He / ONe) |
| 18 | Neutron star |
| 19 | Black hole |

---

## Filtering the output

> The bundled `rdc_SeBa` post-processing tool **does not read this fork's
> 9-column output** — its reader still expects upstream's 18-field format, so it
> silently misaligns every field. Do not use it here. See
> [`README_rdc_SeBa.md`](README_rdc_SeBa.md) for details.

Because `SeBa.data` is 9 whitespace-separated columns with exactly two lines per
binary (T=0 and T=end), filtering is a few lines of Python. Example — final state
of all double white dwarf (types 12–14) systems:

```python
WD = {12, 13, 14}
finals = {}
for line in open("SeBa.data"):
    c = line.split()
    bid, t = c[0], float(c[2])
    if bid not in finals or t >= float(finals[bid][2]):
        finals[bid] = c
for c in finals.values():
    if int(c[5]) in WD and int(c[7]) in WD:   # cols 6 & 8: stellar types
        print(" ".join(c))
```

Adjust the type set and column indices (binary type in col 2, masses in cols 7
and 9) for other selections. This selects on stellar type alone, so it also
catches *disrupted* pairs (binary type 7) that are no longer bound. To keep only
bound systems, additionally require the binary type — e.g. add
`int(c[1]) == 2 and` (2 = Detached) to the `if`.

---

## Performance

This fork replaces the per-event, per-binary I/O of upstream with a compact
two-line-per-binary format and a single output-file handle, and hoists hot-path
math out of `pow()`. Measured single-threaded on macOS against a 67,045-binary
population (`-T 12550 -s 42`):

| | Pre-optimisation (`271c11d`) | Current | Δ |
|-|------------------------------|---------|---|
| User time | 290s | 228s | −21% |
| Sys time  | 66s  | 1s   | −98% |
| Real time | 365s | ~230s | −37% |

The sys-time collapse comes from suppressing terminal spam and opening the output
file once instead of per binary. On the cluster the OpenMP loop scales the
evolution phase roughly as `T_single / N_cores`. Full benchmark details and the
optimisation list are in [`CLAUDE.md`](CLAUDE.md#performance-optimisations).

---

## Known caveats

- **Absolute paths for `-I`**: an unopenable input path silently infinite-loops.
- **Seeded runs need `OMP_NUM_THREADS=1`**: the RNG is not thread-safe, so kick
  angles (and therefore results) are non-deterministic above one thread.
- **`double_star::suppress_output` must stay `true`**: flipping it reactivates the
  old per-binary file writes and defeats the single-file-open optimisation.
- **Longstanding integer-division bug** in `de_dt_gwr()` (`double_star.C`): an
  eccentricity correction term is dropped due to integer division. Present in
  upstream too — see [`CLAUDE.md`](CLAUDE.md#known-issues-and-caveats) before
  touching it, as fixing it invalidates the accuracy-benchmark reference.

---

## References

- Portegies Zwart S.F. & Verbunt F., 1996, A&A, 309, 179:
  *Population synthesis of high-mass binaries*
- Toonen S., Nelemans G., Portegies Zwart S.F., 2012, A&A, 546A, 70T:
  *Supernova Type Ia progenitors from merging double white dwarfs*
</content>
</invoke>
