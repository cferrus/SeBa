# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build

SeBa requires only a C++ compiler. Build in two stages — the top-level step builds the shared libraries, then `dstar/` builds the executable:

```bash
make clean && make   # builds libsstar.a, libnode.a, libstd.a, librdc.a in their subdirs
cd dstar && make     # builds libdstar.a and the SeBa executable
```

The executable lands at `dstar/SeBa`. `Makefile.inc` is auto-generated from `Makefile.inc.conf` by substituting the absolute path; delete it to force regeneration.

**Important build gotcha:** The `dstar/SeBa` executable only relinks when `SeBa.C` itself is newer than the binary. After modifying library code (e.g. `double_star.C`), force a relink manually:
```bash
rm -f dstar/SeBa && cd dstar && make SeBa
```

There are no automated tests. To verify a build, run a quick single-system evolution:

```bash
cd dstar
./SeBa -M 2 -m 1 -e 0.2 -a 200 -T 13500 -z 0.001
```

Output goes to `SeBa.data` in the current directory. Expected: exactly 2 lines (T=0 and T=end).

## Running SeBa

Three modes:

| Mode | Example |
|---|---|
| Single system | `./SeBa -M 2 -m 1 -e 0.2 -a 200 -T 13500 -z 0.001` |
| Input file | `./SeBa -I /path/to/input.txt -T 12550` (format: `a e M m z` per line) |
| Random population | `./SeBa -R -n 250000 -T 13500 -f 4 -m 0.95 -M 10` |

The large benchmark input file lives at `/Users/melz/Work_Program/Data/SeBa_input_T_12550.txt` (67,045 binaries). Always use the full absolute path — relative paths cause a silent infinite loop in `read_binary_params` if the file isn't found.

## Output format

`SeBa.data` is written with **9 columns only** (trimmed from the original 18). Each binary produces exactly 2 lines: T=0 (initial state) and T=end (final state). No intermediate steps are written.

| Col | Field |
|-----|-------|
| 1 | binary identity |
| 2 | binary type (enum) |
| 3 | time (Myr) |
| 4 | semi-major axis (Rsun) |
| 5 | eccentricity |
| 6 | primary stellar type (integer) |
| 7 | primary mass (Msun) |
| 8 | secondary stellar type (integer) |
| 9 | secondary mass (Msun) |

Stellar type integers: 3=Main_Sequence, 5=Hertzsprung_Gap, 6=Sub_Giant, 7=Horizontal_Branch, 8=Super_Giant, 10=Helium_Star, 11=Helium_Giant, 12–14=White_Dwarf variants, 18=Neutron_Star, 19=Black_Hole.

## Output suppression mechanism

All intermediate file I/O is suppressed via a static flag on `double_star`:

- `double_star::suppress_output` — static bool, set to `true` in `main()` unconditionally
- `dump(char*, bool)` — returns immediately if `suppress_output` is true; catches all writes from `sstar/` stellar-type classes (neutron_star, black_hole, etc.) that call `get_binary()->dump()`
- `dump_unconditional(char*, bool)` — bypasses the flag; was used for T=0 and T=end writes before the ostream refactor
- `dump(ostream&, bool)` — the actual formatter; no suppression check (suppression is at the char* entry point)
- `binev.data` is never created (all writes to it are suppressed)

In `main()`, `SeBa.data` is opened **once** as an `ofstream` and passed as `ostream&` to `evolve_binary()`, which calls `ds->dump(outstream, true)` directly. This eliminates ~134K open/close syscalls for a 67K-binary run.

Unguarded `cerr` diagnostic messages (spiral-in, common envelope, Roche contact events) are also gated on `!suppress_output`.

## Performance optimizations applied

Committed on branch `main` of `cferrus/SeBa`:

1. **`cbrt()` for cube roots** — both `roche_radius()` overloads use `cbrt(mr)` instead of `pow(mr, 1/3)` and `q1_3*q1_3` instead of `pow(q1_3, 2)`. Hot path: called on every timestep refinement in `recursive_binary_evolution`.

2. **`pow()` → multiplication** in angular momentum loss functions:
   - `pow(sma, 4)` → `sma2*sma2` in `gwr_angular_momentum_loss()`
   - `pow(semi, 5)` → `semi2*semi2*semi` in `mb_angular_momentum_loss()`
   - `pow(get_total_mass(), 2)` → `m_tot*m_tot`

3. **Cached physical constants** — `c_gwr` and `c_mb` are `static const` locals; computed once from `cnsts` physical constants (which are global and never change) instead of on every timestep call.

4. **Suppressed cerr/cout spam** — ~67K `cout` writes in `add_star.C` commented out; unguarded `cerr` events in `double_star.C` gated on `!suppress_output`.

5. **Single file open/close** — `SeBa.data` opened once in `main()`, eliminating per-binary file I/O overhead.

**Benchmark result** (67,045 binaries, `-T 12550`):
| Metric | Before optimizations | After | Δ |
|--------|---------------------|-------|---|
| Wall time | ~264s | ~238s | -10% |
| User time | ~247s | ~230s | -7% |
| Sys time | ~14.4s | ~7.5s | -48% |

## Accuracy benchmark

Any change to the physics or output logic must be validated by comparing `SeBa.data` against the reference file:

```
/Users/melz/Work_Program/Data/clean12550.data
```

This is the known-good output for the 67,045-binary input run with `-I /Users/melz/Work_Program/Data/SeBa_input_T_12550.txt -T 12550`. It has 134,090 lines (2 per binary). To verify a new build:

```bash
cd dstar
./SeBa -I /Users/melz/Work_Program/Data/SeBa_input_T_12550.txt -T 12550 2>/dev/null
diff SeBa.data /Users/melz/Work_Program/Data/clean12550.data
```

A clean diff (no output) means the new build is numerically identical to the reference. Any differences indicate a regression in physics or output formatting and must be investigated before committing.

## Filtering output with rdc_SeBa

`rdc/rdc_SeBa` filters `SeBa.data` by stellar type and binary state. Must be built as part of the top-level `make`:

```bash
less SeBa.data | ./rdc/rdc_SeBa -f -p white_dwarf -s white_dwarf > SeBa_wdwd.data
```

Key flags: `-f` (first occurrence only), `-R` (full evolution), `-p`/`-s` (primary/secondary type), `-B` (binary type).

## Architecture

The codebase is a C++ class hierarchy rooted at `starbase` → `star`:

- **`star`** (`include/star/star.h`): abstract base for all evolving objects; holds the `node*` tree pointer and unit-conversion factors.
- **`single_star`** (`include/star/single_star.h`): concrete stellar evolution. Each stellar phase is a derived class (`main_sequence`, `hertzsprung_gap`, `sub_giant`, `horizontal_branch`, `super_giant`, `helium_star`, `helium_giant`, `white_dwarf`, `neutron_star`, `black_hole`, etc.) living in `sstar/starclass/`. Phase transitions work by replacing the `single_star` object on the `node`.
- **`double_star`** (`include/star/double_star.h`, `dstar/starclass/double_star.C`): wraps a binary system. Owns the orbital parameters (`semi`, `eccentricity`, `bin_type`) and drives timestep-coupled evolution of both components. Mass-transfer stability, common-envelope events, and tidal circularisation are handled here.
- **`node`** / **`dyn`** (`node/`, `include/`): the N-body tree structure inherited from Starlab. In SeBa only the tree bookkeeping is used; direct gravitational dynamics are not exercised.
- **`cnsts`** (`sstar/starclass/constants.C`): a global `stellar_evolution_constants` object. Holds physical constants (`solar_mass`, `solar_radius`, `G`, `C`, `Myear`, etc.) and simulation parameters (`corotation_eccentricity`, `magnetic_braking_exponent`, etc.). Values are fixed at startup and never change per-binary.

### Directory map

```
include/          Global headers (star/, node-level)
sstar/starclass/  Single-star phase implementations (one .C per stellar type)
dstar/starclass/  double_star and binary support
dstar/stardyn/    Binary-state evolution logic (mass transfer, CE)
dstar/init/       Initialisation of double systems
dstar/util/       dstar-level utilities (dump, semi↔period conversion)
dstar/io/         Binary I/O
dstar/SeBa.C      Main entry point; CLI parsing, population loop
rdc/              Post-processing tools (rdc_SeBa and friends)
std/              Shared low-level utilities
node/dyn/         Kepler orbit routines used by double_star
starev.C          Standalone single-star evolution tool (top-level)
```

### Key compile-time flags

- `-DTOOLBOX`: activates the `main()` in each `.C` file (allows per-file standalone builds).
- `-DHAVE_CONFIG_H`: enables `config.h` inclusion.
- Debug verbosity is controlled by `#define REPORT_*` booleans near the top of `double_star.C` (`REPORT_BINARY_EVOLUTION`, `REPORT_RECURSIVE_EVOLUTION`, `REPORT_FUNCTION_NAMES`, `REPORT_TRANFER_STABILITY`) — all `false` by default; flip to `true` for detailed trace output.

### Pending work

- **OpenMP parallelization** of the binary loop in `dstar/SeBa.C` — each binary is independent so the loop is embarrassingly parallel. Expected near-linear speedup with core count.
