# CLAUDE.md

This file documents the SeBa codebase as it lives in the `cferrus/SeBa` fork. It is the primary reference for build instructions, output format, performance decisions, accuracy validation, and cluster usage. Keep it up to date when making changes.

---

## Table of contents

1. [Build](#build)
2. [Running SeBa](#running-seba)
3. [Cluster workflow](#cluster-workflow)
4. [Output format](#output-format)
5. [Output suppression mechanism](#output-suppression-mechanism)
6. [OpenMP parallelisation](#openmp-parallelisation)
7. [Performance optimisations](#performance-optimisations)
8. [Accuracy benchmark](#accuracy-benchmark)
9. [Filtering output with rdc_SeBa](#filtering-output-with-rdc_seba)
10. [Architecture](#architecture)
11. [Known issues and caveats](#known-issues-and-caveats)

---

## Build

SeBa requires only a C++ compiler and `make`. No external dependencies.

### Standard build (macOS / local)

```bash
# From the repo root:
make clean && make        # builds libsstar.a, libnode.a, libstd.a, librdc.a
cd dstar && make          # builds libdstar.a and the SeBa executable
```

The executable lands at `dstar/SeBa`.

On macOS with Apple Clang, OpenMP is not available without `libomp`. The Makefile detects this automatically — the binary compiles and runs correctly in single-threaded mode.

### Cluster build (Linux / GCC)

Identical commands. On Linux the `dstar/Makefile` detects `uname -s == Linux` and adds `-fopenmp` to both compile and link steps, giving full multi-core parallelism automatically.

```bash
make clean && make
cd dstar && make
```

To verify OpenMP is active after building on the cluster:

```bash
ldd dstar/SeBa | grep omp     # should show libgomp or similar
```

### Makefile.inc regeneration

`Makefile.inc` is auto-generated from `Makefile.inc.conf` by substituting the absolute repo path. It is not committed. Delete it to force regeneration:

```bash
rm -f Makefile.inc dstar/Makefile.inc && make clean && make
```

### Force-relink after library changes

The `dstar/SeBa` executable only relinks when `SeBa.C` itself is newer than the binary. After editing library code (e.g. `double_star.C`, `single_star.C`), force a full relink:

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

### Modes

| Mode | Command |
|------|---------|
| Single system | `./SeBa -M 2 -m 1 -e 0.2 -a 200 -T 13500 -z 0.001` |
| Input file | `./SeBa -I /path/to/input.txt -T 12550` |
| Random population | `./SeBa -R -n 250000 -T 13500 -f 4 -m 0.95 -M 10` |

Input file format: one binary per line as `a  e  M  m  z` (semi-major axis, eccentricity, primary mass, secondary mass, metallicity).

**Always use absolute paths for `-I`.** A relative path that cannot be opened causes a silent infinite loop in `read_binary_params`.

### Key flags

| Flag | Meaning |
|------|---------|
| `-I file` | Read initial conditions from file |
| `-R` | Generate random initial conditions |
| `-n N` | Number of binaries (random mode) |
| `-N offset` | Starting binary identity number |
| `-T Myr` | End time in Myr |
| `-s seed` | Fix random seed (required for reproducibility) |
| `-z Z` | Metallicity (0.0001–0.03) |
| `-O file` | Output filename (default: `SeBa.data`) |
| `-V` | Verbose: print diagnostics to stderr/stdout |

### Verbose flag

By default SeBa produces **no terminal output**. Add `-V` to enable diagnostics:

```bash
./SeBa -V -M 2 -m 1 -e 0.2 -a 200 -T 13500 -z 0.001
```

With `-V`:
- **stdout**: run log story (seed, CE/kick settings, distribution functions) printed at end of `main()`
- **stderr**: per-component phase transitions, merger/spiral-in events, kick velocity diagnostics, physics warnings

Omit `-V` for all batch and pipeline runs.

---

## Cluster workflow

The cluster run is driven by `Data/run.sh`, which calls `SeBaCluster` (the SeBa binary, renamed or symlinked) once per epoch/population slice. There are ~13,500 jobs, each writing to a unique output file.

### Prerequisites

- `SeBaCluster` binary must be present in the `Data/` directory (copy or symlink `dstar/SeBa`)
- Python 3 with `psutil` installed (`pip install psutil`)
- Input files (`SeBa_input_T_*.txt`) must be present in `Data/`

### Running in parallel

```bash
cd Data/
python3 run_seba.py
```

`run_seba.py` handles everything:
1. Reads `run.sh` and rewrites each job's `-O trial.data` to a unique `-O trial_{i}.data`
2. Runs all jobs in parallel using `ThreadPoolExecutor`, with worker count dynamically adjusted every 3 minutes based on live CPU load (using `psutil`)
3. Merges all `trial_*.data` files into `trial.data` and removes the per-job files

Worker count is bounded between `MIN_WORKERS=4` and `MAX_WORKERS=54`; edit those constants at the top of the script to match your machine.

If `psutil` is unavailable, `run_parallel.sh` provides a GNU `parallel` fallback:

```bash
./run_parallel.sh 32    # replace 32 with your core count
```

### Why not the original run.sh directly?

`run.sh` passes `-O trial.data` to every job. `SeBa.C` opens that file with `ios::trunc`, which **truncates and overwrites** on each invocation. Running jobs in parallel would cause every process to destroy each other's output. Running sequentially would leave only the last job's results. Both `run_seba.py` and `run_parallel.sh` fix this by assigning unique output filenames and merging at the end.

### Copying the binary to Data/

```bash
cp dstar/SeBa Data/SeBaCluster
# or:
ln -sf ../SeBa/dstar/SeBa Data/SeBaCluster
```

---

## Output format

`SeBa.data` is written with **9 columns** (trimmed from the original 18-column format). Each binary produces exactly **2 lines**: T=0 (initial state) and T=end (final state). No intermediate timesteps are written.

| Col | Field | Units |
|-----|-------|-------|
| 1 | binary identity | integer |
| 2 | binary type (enum) | integer |
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
| 3 | Semi-detached |
| 7 | Disrupted |
| 9 | Merged |

### Stellar type integers

| Value | Type |
|-------|------|
| 3 | Main_Sequence |
| 5 | Hertzsprung_Gap |
| 6 | Sub_Giant |
| 7 | Horizontal_Branch |
| 8 | Super_Giant |
| 10 | Helium_Star |
| 11 | Helium_Giant |
| 12–14 | White_Dwarf variants |
| 18 | Neutron_Star |
| 19 | Black_Hole |

---

## Output suppression mechanism

Three layers work together to keep SeBa silent by default.

### `double_star::suppress_output`

- Static bool, set to `true` in `main()` **unconditionally** (not toggled by `-V`).
- `dump(char*, bool)` returns immediately when true — suppresses all `sstar/` stellar-type classes (neutron_star, black_hole, etc.) that call `get_binary()->dump()`.
- `dump(ostream&, bool)` — the actual formatter, no suppression check. This is the path used for deliberate T=0/T=end writes.
- `binev.data` is never created.
- Merger/spiral-in `cerr` diagnostics in `double_star.C` are gated on `!suppress_output`.

**Critical invariant:** `double_star::suppress_output` must remain `true`. If set to `false`, the old per-binary `dump("SeBa.data", true)` call sites would reactivate, duplicating output and defeating the single-file-open optimisation.

### `single_star::suppress_output`

- Static bool in `sstar/starclass/single_star.C`, declared in `include/star/single_star.h`.
- Set to `!verbose` in `main()`.
- Gates phase-transition `cerr` in `star_transformation_story()` and core-mass warnings.

### `cerr` redirect in `main()`

```cpp
static ofstream null_stream("/dev/null");
if (!verbose) cerr.rdbuf(null_stream.rdbuf());
```

Catch-all for any remaining unguarded `cerr` calls anywhere in the codebase (kick diagnostics in `constants.C`, warnings in `neutron_star.C`, `black_hole.C`, etc.).

### Summary table

| Mechanism | Default | With `-V` | Controls |
|-----------|---------|-----------|----------|
| `double_star::suppress_output = true` | active | active | File writes via `dump(char*,bool)` + merger/CE `cerr` |
| `single_star::suppress_output = true` | active | inactive | Phase-transition `cerr` + core-mass warnings |
| `cerr → /dev/null` | active | inactive | All remaining unguarded `cerr` |
| `log_root->print_log_story(cout)` | suppressed | shown | Run log on stdout |

---

## OpenMP parallelisation

### How it works

`dstar/SeBa.C` parallelises the binary evolution loop with OpenMP. The key design decisions:

1. **Pre-read phase**: all binary parameters from the input file are loaded into a `std::vector<BinaryInput>` before the parallel section. This removes the serial `ifstream` bottleneck and allows `n` to be known upfront.

2. **Thread-local trees**: every thread creates its own `dyn` tree (`root_i`, `the_binary_i`) and `double_star` object. No tree state is shared between iterations.

3. **Thread-local output buffering**: each thread writes its T=0 and T=end rows into a `std::ostringstream`. The completed buffer is flushed to the shared `ofstream` inside `#pragma omp critical(output)`. This keeps both rows for a given binary contiguous and avoids interleaved output.

4. **RNG serialisation for `-R` mode**: `mkrandom_binary` uses a global RNG. In random-population mode, the random draw is wrapped in `#pragma omp critical(rng)`. Evolution (the slow part) still runs in parallel; only the draw is serialised.

5. **`schedule(dynamic, 1)`**: binary evolution time varies widely (simple detached pairs vs. common-envelope systems). Dynamic scheduling with chunk size 1 maximises load balance.

### Platform behaviour

| Platform | OpenMP active? | How |
|----------|---------------|-----|
| Linux / GCC | Yes | `dstar/Makefile` detects `uname -s == Linux`, adds `-fopenmp` |
| macOS / Apple Clang | No | Falls back to single-threaded via `#ifdef _OPENMP` guards |
| macOS / Homebrew GCC or LLVM | Yes | Override: `CXX=g++-14 make` |

### Thread safety scope

The following are safe under the current implementation:
- `cnsts` (global `stellar_evolution_constants`): read-only after `main()` initialises it — safe.
- `double_star::suppress_output`, `single_star::suppress_output`: set once before the loop, never written during evolution — safe.
- `cerr` redirect to `/dev/null`: set once before the loop — safe.
- Heap allocation (`mkdyn`, `new_double_star`, etc.): C++ `operator new` is thread-safe — safe.

The following are **not** thread-safe and are handled:
- Global RNG (`srandinter`/`randinter`): serialised via `#pragma omp critical(rng)` in `-R` mode. In `-I` (input file) mode, RNG is only called for NS/BH kick angles — this introduces non-determinism under OpenMP, but `run.sh` does not use `-s`, so results are already non-deterministic between runs. The `-s 42` accuracy benchmark must be run single-threaded (`OMP_NUM_THREADS=1`) for reproducibility.

---

## Performance optimisations

All optimisations are committed on `main` of `cferrus/SeBa`. Listed in chronological order.

### 1. `cbrt()` for cube roots (`dstar/starclass/double_star.C`)

Both `roche_radius()` overloads use `cbrt(mr)` instead of `pow(mr, 1.0/3.0)`, and `q1_3*q1_3` instead of `pow(q1_3, 2)`. Called on every timestep refinement inside `recursive_binary_evolution` — the hottest path in the codebase.

### 2. `pow()` → multiplication in angular momentum loss (`double_star.C`)

- `pow(sma, 4)` → `sma2*sma2` in `gwr_angular_momentum_loss()`
- `pow(semi, 5)` → `semi2*semi2*semi` in `mb_angular_momentum_loss()`
- `pow(get_total_mass(), 2)` → `m_tot*m_tot`

### 3. Cached physical constants (`double_star.C`)

`c_gwr`, `c_mb`, `G3M3_C5R4`, and `G3M3_C5` are all `static const` locals initialised once from the global `cnsts` object. Previously recomputed from `pow()` on every call. These functions are called on every timestep for any binary undergoing angular momentum loss or gravitational-wave inspiral.

### 4. Suppressed cerr/cout I/O spam

- ~67K `cout` writes in `sstar/init/add_star.C` commented out
- Merger/spiral-in `cerr` in `double_star.C` gated on `!suppress_output`
- `PRC(pk)` / `PRL(v_disp)` debug macros in `constants.C` commented out
- All remaining unguarded `cerr` silenced by redirecting to `/dev/null` in `main()`

### 5. Single file open/close (`dstar/SeBa.C`)

`SeBa.data` (or the `-O` output) is opened once as an `ofstream` and passed by reference to `evolve_binary()`. Previously, each binary opened and closed the file for every write, producing ~134K syscalls for a 67K-binary run.

### 6. `pow()` → multiply in `gravrad()` / `de_dt_gwr()` (`double_star.C`)

- `pow(semi, 4.)` → `semi2*semi2`
- `pow(1-e*e, 2.5)` → `beta2*sqrt(beta)` where `beta = 1-e*e`
- `pow(eccentricity, 2.)` → precomputed `ecc2` (appeared 4× in `gravrad()`)

These are on the GWR inspiral code path, active for tight binaries.

### 7. OpenMP parallel binary loop (`dstar/SeBa.C`, `dstar/Makefile`)

See [OpenMP parallelisation](#openmp-parallelisation) above. On the cluster, wall time scales as approximately `T_single / N_cores` for the evolution loop.

### Benchmark results

All measured single-threaded on macOS, 67,045 binaries, `-I SeBa_input_T_12550.txt -T 12550 -s 42`:

| | Original | After 1–5 | After 6–7 | Δ total |
|-|----------|-----------|-----------|---------|
| User time | ~264s | ~230s | ~228s | −14% |
| Sys time  | ~14.4s | ~7.5s | ~1.0s | −93% |

The sys-time drop is dominated by eliminating per-binary file I/O (opt 5) and the `ostringstream` buffering (opt 7).

---

## Accuracy benchmark

Any change to physics or output logic must be validated against the reference file before committing.

### Reference file

```
/Users/melz/Work_Program/Data/12550_cleaned.data
```

- 18-column original SeBa format (not our 9-column format)
- 134,090 lines (2 per binary, 67,045 binaries)
- Produced with commit `af8654c`, `-I SeBa_input_T_12550.txt -T 12550 -s 42`

### Running the check

```bash
cd dstar
OMP_NUM_THREADS=1 ./SeBa -I /Users/melz/Work_Program/Data/SeBa_input_T_12550.txt -T 12550 -s 42
python3 -c "
our, ref = {}, {}
for line in open('SeBa.data'):
    c = line.split(); bid, t = c[0], float(c[2])
    our[(bid, 'T0' if t==0 else 'Tend')] = c
for line in open('/Users/melz/Work_Program/Data/12550_cleaned.data'):
    c = line.split(); bid, t = c[0], float(c[3])
    ref[(bid, 'T0' if t==0 else 'Tend')] = c
mm = [(k[0], o, ref[k]) for k, o in our.items()
      if k[1]=='Tend' and k in ref and
      any([o[1]!=ref[k][1], o[3]!=ref[k][4], o[4]!=ref[k][5],
           o[5]!=ref[k][7], o[6]!=ref[k][8], o[7]!=ref[k][13], o[8]!=ref[k][14]])]
print(f'Mismatches: {len(mm)} / {len(our)//2}')
if mm:
    for bid, o, r in mm[:5]:
        print(f'  {bid}: ours={o[1:][:8]}')
        print(f'       ref ={[r[1],r[4],r[5],r[7],r[8],r[13],r[14]]}')
"
```

**`OMP_NUM_THREADS=1` is required** — the RNG is not thread-safe and produces different kick angles under OpenMP, making the count non-deterministic above 1 thread.

### Expected result

**88 mismatches / 67,045 binaries (0.13%).**

These are a known, stable divergence — small eccentricity differences in WD+Brown Dwarf and WD+companion systems caused by minor I/O code-path differences that shift the global RNG sequence. No `bin_type`, `s1_type`, or `s2_type` differences are expected. A count above 88, or any new type differences, indicates a regression.

---

## Filtering output with rdc_SeBa

`rdc/rdc_SeBa` filters `SeBa.data` by stellar type and binary state. Built as part of the top-level `make`.

```bash
# WD+WD systems only:
less SeBa.data | ./rdc/rdc_SeBa -f -p white_dwarf -s white_dwarf > SeBa_wdwd.data
```

Key flags:

| Flag | Meaning |
|------|---------|
| `-f` | First occurrence only (initial state snapshot) |
| `-R` | Full evolution history |
| `-p type` | Filter on primary stellar type |
| `-s type` | Filter on secondary stellar type |
| `-B type` | Filter on binary type |

---

## Architecture

The codebase is a C++ class hierarchy rooted at `starbase` → `star`, inherited from the Starlab N-body framework. Only the stellar evolution and binary interaction components are actively used; the gravitational dynamics engine is not exercised.

### Core classes

- **`star`** (`include/star/star.h`): abstract base for all evolving objects. Holds the `node*` tree pointer and unit-conversion factors.

- **`single_star`** (`include/star/single_star.h`): concrete stellar evolution. Each stellar phase is a separate derived class (`main_sequence`, `hertzsprung_gap`, `sub_giant`, `horizontal_branch`, `super_giant`, `helium_star`, `helium_giant`, `white_dwarf`, `neutron_star`, `black_hole`, etc.) in `sstar/starclass/`. Phase transitions replace the `single_star` object on the node.

- **`double_star`** (`include/star/double_star.h`, `dstar/starclass/double_star.C`): wraps a binary system. Owns orbital parameters (`semi`, `eccentricity`, `bin_type`) and drives timestep-coupled evolution of both components. Mass-transfer stability, common-envelope events, tidal circularisation, gravitational-wave inspiral, and magnetic braking are handled here.

- **`node`** / **`dyn`** (`node/`, `include/`): N-body tree structure from Starlab. Used only for tree bookkeeping in SeBa.

- **`cnsts`** (`sstar/starclass/constants.C`): global `stellar_evolution_constants` object. Holds all physical constants (`solar_mass`, `solar_radius`, `G`, `C`, `Myear`, etc.) and simulation parameters (`corotation_eccentricity`, `magnetic_braking_exponent`, CE efficiency, etc.). Initialised once at startup, never modified per-binary — read-only and thread-safe.

### Hot code paths

In order of CPU time consumed (approximate, unprofilesd):

1. `double_star::recursive_binary_evolution()` — adaptive timestepping loop; calls `evolve_element()` and `roche_radius()` on every step.
2. `double_star::roche_radius()` — called multiple times per step; uses `cbrt()` (optimised).
3. `double_star::gwr_angular_momentum_loss()` / `mb_angular_momentum_loss()` — called every step for interacting binaries.
4. `single_star` phase evolution functions — one per stellar type, called via virtual dispatch.

### Directory map

```
include/           Global headers (star/, node-level)
sstar/starclass/   Single-star phase implementations (one .C per stellar type)
dstar/starclass/   double_star and binary orbital mechanics
dstar/stardyn/     Binary-state evolution (mass transfer, common envelope)
dstar/init/        Initialisation of double systems
dstar/util/        dstar-level utilities (dump, semi↔period conversion)
dstar/io/          Binary I/O
dstar/SeBa.C       Main entry point: CLI parsing, pre-read, parallel binary loop
rdc/               Post-processing tools (rdc_SeBa and friends)
std/               Shared low-level utilities (RNG, I/O helpers)
node/dyn/          Kepler orbit routines used by double_star
starev.C           Standalone single-star evolution tool (top-level)
Data/run.sh        Sequential job list (~13,500 SeBaCluster invocations)
Data/run_parallel.sh  Parallel wrapper for run.sh using GNU parallel
Data/slurm_array.sh   SLURM job array template (if applicable)
```

### Key compile-time flags

| Flag | Effect |
|------|--------|
| `-DTOOLBOX` | Activates `main()` in each `.C` file for standalone per-file builds |
| `-DHAVE_CONFIG_H` | Enables `config.h` inclusion |
| `_OPENMP` | Defined by compiler when `-fopenmp` is active; guards OpenMP code in `SeBa.C` |

Debug verbosity: flip the `#define REPORT_*` booleans near the top of `double_star.C` (`REPORT_BINARY_EVOLUTION`, `REPORT_RECURSIVE_EVOLUTION`, `REPORT_FUNCTION_NAMES`, `REPORT_TRANFER_STABILITY`) from `false` to `true` for a detailed trace.

---

## Known issues and caveats

### Integer division in `de_dt_gwr()` (`double_star.C:2052`)

```cpp
real de_dt = (-305/15)*e*(1 + (121/304)*e*e)/...
```

`(-305/15)` evaluates to `-20` (integer division; correct value is `-20.33...`) and `(121/304)` evaluates to `0` (correct value is `~0.398`), completely dropping the eccentricity correction term. This is a longstanding bug present in the upstream code and in the `af8654c` reference. **Do not fix it without regenerating the reference file**, as doing so would shift the mismatch baseline and invalidate all prior accuracy comparisons.

### RNG is not thread-safe under OpenMP

The global RNG used for NS/BH kick angles produces data races when the binary loop runs with `OMP_NUM_THREADS > 1`. This introduces non-determinism in kick velocities. For `-I` (input file) mode this is acceptable — `run.sh` does not use `-s`, so results are already non-deterministic between runs. For reproducible single-binary or seeded runs, always set `OMP_NUM_THREADS=1`.

### `run.sh` truncation bug

Every job in `run.sh` passes `-O trial.data`. `SeBa.C` opens this file with `ios::trunc`, meaning each sequential invocation overwrites the previous output. Use `run_parallel.sh` instead, which assigns unique output filenames and merges at the end.

### `dstar/SeBa` relink gotcha

After editing any library file (anything outside `dstar/SeBa.C`), the executable does not automatically relink. Run:
```bash
rm -f dstar/SeBa && cd dstar && make SeBa
```

### Profiling next steps

The 228s single-threaded wall time has not been formally profiled. The bottleneck almost certainly lives in `recursive_binary_evolution` → `evolve_element` → `roche_radius`. To profile:

```bash
# Rebuild with profiling:
cd dstar/starclass && touch double_star.C && make CXXFLAGS+=-pg
cd .. && make CXXFLAGS+=-pg SeBa
# Run:
./SeBa -I /path/to/input.txt -T 12550 -s 42
gprof SeBa gmon.out | head -40
```
