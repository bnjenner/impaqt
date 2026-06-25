# C++ Cleanup — Roadmap & Status

Working branch: **`cpp-audit-fixes`** (off `master`, started 2026-06-25). Nothing
pushed or merged yet. 20 commits so far, each built clean under `-Wall -Wextra`
(C++17) and verified before committing.

Goal: incremental C++ cleanup from a codebase audit. One logical change per commit;
behavior preserved unless a change is explicitly a bug fix.

---

## ✅ Done (on the branch)

**Correctness bugs**
- Uninitialized `bool` in `ContainmentList::collapse_intervals` (UB) — `0d87a83`
- `ClusterList::get_transcript_num` accumulated into a member (double-count) — `bcba97c`
- `thread_queue::dispatch(const call&&)` + `std::move` silently copied — `5b22657`
- `throw "string literal"` → `std::runtime_error` + try/catch in `main` — `1f57c4b`
- **`paths` string-map → `struct Path` — fixed a latent bug**: cluster indices ≥10
  were mangled by the 2-char-string encoding. Found via the real-data guard at
  `chr1:160862023-160866444` (n5=12, n3=10), the only such locus in 3.8M reads — `2b5cae6`

**Performance**
- Heavy getters returned vectors by value → return by `const&` (killed O(n²)
  per-iteration copies in the assignment loop) — `ab916b0`

**Memory / RAII**
- `processes` now `vector<unique_ptr<Impaqt>>` (also fixes a leak on no-annotation runs) — `6e31cf9`

**Build / warnings**
- C++14 → C++17, `-Wall -Wextra` on first-party targets only — `33138b5`
- `-Wsign-compare` cleanup (signed indices / cache size as `const int`) — `bd1016e`
- Dead locals (`243e2a4`) and dead parameters (`16e3349`) from `-Wextra`

**Idioms / hygiene**
- `NULL` → `nullptr` — `b7f9cd0`
- `variable_swap` → `std::swap`; one-line `file_exists` — `f43abe1`
- Deduped `GeneNode::overlap` into `check_bounds` — `520deef`
- `std::endl` → `"\n"` for end-of-output — `e9f5376`
- Header hygiene: direct `global_args.h` include, `inline argparse`, `set()` vs `file(GLOB)` — `671735b`
- `const`-qualified read-only getters on all four node/list classes — `34623d2`, `1368fdf`

**Repo**
- Removed accidentally-nested duplicate `test/test/` trees — `65b0d3d`
- Ignore `temp/` (local benchmark fixtures) — `b0ec403`

---

## ⏳ Remaining from the audit

| Item | Notes | Effort | Risk |
|---|---|---|---|
| **RAII for the linked lists** | `ClusterList`/`AnnotationList` still raw `new`/`delete` with hand-written teardown. Owning `unique_ptr` nodes, or a flat pool. Remove the no-op `~ClusterNode`/`~GeneNode`. | Medium | Watch destructor recursion depth on long lists. |
| **Manual 1000-elem resize** | `ClusterNode::add_alignment` reimplements vector growth with a `% 1000` branch → `reserve` + `push_back`. | Low | Low |
| **Deeper const-threading** | Getters are `const` now; read-only free funcs (`get_read_overlap`, `get_transcript_overlap`, `assign_*`, etc.) still take non-`const` `GeneNode*`/`ClusterNode*`. Thread `const` through. | Medium | Low, but ripples through signatures. |
| **Assignment unit test** | `assign_*` logic has no dedicated unit test (the old stub was an empty `TEST_F`). | Medium | — |
| **≥10-cluster regression test** | Pin the `paths` bug we fixed: a `get_linked_clusters`/`get_coordinates` test with cluster indices ≥10. | Low–Med | — |

---

## 📋 Separately tracked (raised during the session)

### Dependency stack (large)
Replace the vendored `ext/` copies (~1200 committed files).
- **Drop seqan entirely.** Used only in `include/ArgParser.h` for `seqan::ArgumentParser`
  + a file-extension helper. Replace the arg parser (hand-rolled, or CLI11 / cxxopts).
  Must re-implement: positional BAM arg; options `-t -a -s -n -q -w -m -p -e -d -f -u -i -o`
  with their current defaults; `.bam`/`.gtf`/`.gff` extension checks.
- **bamtools**: integrate via FetchContent or submodule; likely **update to a recent release**.
- **googletest**: FetchContent/submodule.
- Caveat: deps must still compile under C++17 + `-Wall -Wextra`. The Apple Silicon CI
  job below would surface libc++ friction.

### Apple Silicon CI (small)
Add a `macos-14` (Apple Silicon, arm64) matrix job to `.github/workflows/c-cpp.yml`.
Proposed:
```yaml
jobs:
  build_and_test:
    strategy:
      fail-fast: false
      matrix:
        os: [ubuntu-latest, macos-14]   # macos-14 = Apple Silicon (arm64)
    runs-on: ${{ matrix.os }}
    steps:
    - uses: actions/checkout@v4
    - name: Install deps (macOS)
      if: runner.os == 'macOS'
      run: brew install cmake zlib
    - run: cmake -DCMAKE_BUILD_TYPE=Release -S . -B build
    - run: cmake --build build
    - run: ctest --test-dir build --output-on-failure
```
**Caveat — `long double` width:** x86-64 = 80-bit, Apple Silicon = 64-bit. Expression/
count accumulation uses `long double`, so per-gene counts can differ in the last digits
across architectures (and the md5 guard below is therefore arch-specific). Not a bug;
pin to `double` only if exact cross-platform reproducibility is required.

---

## 🧪 Verification setup (recreate tomorrow)

Local benchmark/regression fixtures live in `temp/` (gitignored, ~1.1 GB, **not** in the repo):
- `temp/TS25_2_1_dedup.bam` (+ `.bai`) — real mouse data, 3.8M reads, 61 contigs, coordinate-sorted, **forward-stranded**.
- `temp/gencode.vM35.annotation.gtf` — GENCODE vM35 (matches the BAM's `chr` naming).
  Re-fetch: `wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M35/gencode.vM35.annotation.gtf.gz` then `gunzip`.

Output is deterministic across reruns **and** thread counts. The guard hashes the GTF
**body only** (excludes `##` header lines, which embed the input path).

Two scripts (were in `/tmp`; copy back if the box rebooted):

### `verify_impaqt.sh` — fast guard (unit tests + small e2e)
Rebuilds, runs `ctest`, and diffs end-to-end output for the small `test/data` BAMs
against a captured baseline. Run after every change.

### `bench_impaqt.sh` — real-data guard + timing
```bash
#!/bin/bash
set -e
ROOT=/home/bnjenner/bnjenner_software/impaqt
BAM=$ROOT/temp/TS25_2_1_dedup.bam
GTF=$ROOT/temp/gencode.vM35.annotation.gtf
THREADS=${1:-16}

# Post-Path-refactor baselines (body-only md5)
BASE_NOANNO_GTF=67db236a2ceef4672ff8bff5ac287046   # transcript identification / paths
BASE_ANNO_GTF=d5eb90128d1fc851893b23e02a340664     # assignment GTF
BASE_ANNO_COUNTS=9a31c65b823dd40c3da94c1ceddb684b  # per-gene counts + summary

fail=0
check() { local got; got=$(grep -v '^##' "$2" | md5sum | cut -d' ' -f1)
  if [ "$got" = "$3" ]; then echo "  OK   $1 ($got)"; else echo "  FAIL $1 (got $got want $3)"; fail=1; fi; }

echo "== transcript identification (no annotation) =="
"$ROOT/build/impaqt" "$BAM" -t "$THREADS" -o /tmp/g_noanno.gtf >/dev/null 2>/dev/null
check "noanno gtf" /tmp/g_noanno.gtf "$BASE_NOANNO_GTF"

echo "== gene assignment (forward + annotation) =="
"$ROOT/build/impaqt" "$BAM" -a "$GTF" -t "$THREADS" -o /tmp/g_anno.gtf >/tmp/g_anno.counts 2>/dev/null
check "anno gtf"    /tmp/g_anno.gtf    "$BASE_ANNO_GTF"
check "anno counts" /tmp/g_anno.counts "$BASE_ANNO_COUNTS"
[ $fail -eq 0 ] && echo "ALL GUARDS PASS" || { echo "REGRESSION"; exit 1; }
```

Reference performance (16 threads): ~27s no-annotation, ~29s with assignment, ~1.4 GB peak RSS.

Assignment sanity (forward, vM35): 2,153,052 assigned / 196,777 unassigned / 10,986
ambiguous / 1,450,774 multimapping (NH>1, excluded by default) / 15,855 transcripts.

---

## 🗺️ Suggested order for tomorrow

1. **Quick wins first**: `ClusterNode::add_alignment` resize → `reserve`; then the
   ≥10-cluster regression test (cheap, and it locks in the bug we just fixed).
2. **Assignment unit test** — gives real coverage before deeper refactors.
3. **RAII for the lists** — do it against both guards; watch destructor recursion.
4. **Deeper const-threading** — mechanical once #3 is in.
5. Then pick up the **dependency stack** (start by dropping seqan / replacing the arg
   parser, since it's self-contained) and add the **Apple Silicon CI** job.

Always: `verify_impaqt.sh` after every change; `bench_impaqt.sh` for anything touching
DBSCAN/paths or the assignment path.
