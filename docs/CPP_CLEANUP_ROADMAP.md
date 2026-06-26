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

**Quick wins (2026-06-26 session)**
- `ClusterNode::add_alignment` manual `% 1000` resize → `reserve` + `push_back` — `dd10a36`
- `≥10`-cluster path regression test (`LargeClusterIndexPaths`, pins `2b5cae6`) — `9d8fc84`
- Gene-assignment unit tests (new `assign_test` target) — `aa9d189`

**Correctness / memory (2026-06-26 session)**
- `ClusterList::delete_list()` null-head crash — guarded; `DestroyEmptyList` regression test (SEGFAULTs without the fix) — `edf0862`
- **Memory-leak audit (ASan/LSan).** Found one real leak: `Impaqt::cluster_list` was a
  raw `ClusterList*`, `new`'d in `create_clusters()`, never freed (empty `~Impaqt`) →
  every contig's list + nodes + transcripts leaked at exit (LSan: single root, 26
  allocs). Fixed by making it `std::unique_ptr<ClusterList>` — `f5640e7`. The rest of
  the raw `new`/`delete` (`merge_nodes`, `delete_list`, `AnnotationList`) is **ASan-clean
  on the real 3.8M-read data** — no leaks, use-after-free, or double-free. This is the
  evidence backing the decision to leave the linked-list internals as raw pointers.

---

## ⏳ Remaining from the audit

| Item | Notes | Effort | Risk |
|---|---|---|---|
| **RAII for the linked lists** | **Decided: leave as-is.** `ClusterList`/`AnnotationList` keep raw `new`/`delete`. The ASan audit proved the teardown + `merge_nodes` are leak/UAF/double-free-clean on real data; a `unique_ptr<next>` rewrite trades that for recursion + merge-ownership churn (see the `merge_with_next` discussion). Revisit only if the lists ever outlive a single process. | — | — |
| **Deeper const-threading** | ✅ **done** `b1e9311`. `get_read_overlap`/`get_transcript_overlap` → `const GeneNode*`; `get_closest_gene` clust + `assign_reads_to_genes` node → `const ClusterNode*`. Mutating funcs correctly left non-const. Covered by `assign_test`. | — | — |

---

## 📋 Separately tracked (raised during the session)

### Dependency stack ✅ DONE
Vendored `ext/` (1201 files, ~30 MB) fully replaced; **ext/ removed entirely.**
- ✅ **Drop seqan** — `5efd998`. Was used only for CLI parsing; replaced `ArgParser.h`
  with a hand-rolled parser (no new dep), `argparse` returns a `ParseStatus` enum.
- ✅ **bamtools → FetchContent v2.5.3** — `774a78f`, pinned by immutable SHA `7e5c183`.
  The vendored copy was a hand-patched ~2.5.1 (std::-qualified C funcs to build on
  modern libstdc++); v2.5.3 already carries that fix upstream, so it builds clean under
  C++17. Full `MakeAvailable` (the `BamTools` lib target lives in `src/` in 2.5.3);
  headers `SYSTEM`. Side effect: 2 bamtools toolkit tests register in ctest and pass.
- ✅ **googletest → FetchContent release-1.12.1** — `5379f43`, pinned by SHA `58d77fa`
  (same version as vendored → pure structural change). Replaced `cmake/gtest.cmake`
  (ExternalProject) + removed `GTEST_INCLUDE_DIRS`/`GTEST_LIBS_DIR` wiring.
- ✅ `ext/` removed — `<sha>` (see git log). Provenance now lives in CMakeLists comments.

**How "correct version" is guaranteed:** each dep pinned to a full commit SHA (a bare
tag can be force-moved; a SHA is content-addressed → reproducible), and verified by
byte-identical guard output (proves identical BAM decoding) + all tests. A fresh
configure re-fetches both in ~9 s.

Caveat for the Apple Silicon job below: deps build under C++17 + `-Wall -Wextra`
(bamtools headers are `SYSTEM`, so its warnings don't fail our build); arm64 libc++
is the remaining unknown.

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

### ASan / LeakSanitizer (memory audit)
No valgrind on this box, but gcc's ASan+LSan works and is faster. Separate build dir:
```bash
cmake -S . -B build-asan -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_CXX_FLAGS="-fsanitize=address -fno-omit-frame-pointer -g" \
  -DCMAKE_EXE_LINKER_FLAGS="-fsanitize=address"
cmake --build build-asan --target impaqt -j4
ASAN_OPTIONS=detect_leaks=1 ./build-asan/impaqt temp/TS25_2_1_dedup.bam \
  -a temp/gencode.vM35.annotation.gtf -t 8 -o /tmp/asan_real.gtf
```
Exit 0 with no `SUMMARY: AddressSanitizer` line = clean. Run the real BAM (not just
`test/data/`) — only it exercises `merge_nodes` and large lists. As of `f5640e7`: clean.

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

## 🗺️ Suggested order

1. ~~Quick wins: `add_alignment` resize → `reserve`; ≥10-cluster regression test.~~ ✅ done `dd10a36`, `9d8fc84`
2. ~~Assignment unit test — real coverage before deeper refactors.~~ ✅ done `aa9d189`
3. **RAII for the lists** — do it against both guards; watch destructor recursion.
4. **Deeper const-threading** — mechanical, now backed by `assign_test`.
5. Then pick up the **dependency stack** (start by dropping seqan / replacing the arg
   parser, since it's self-contained) and add the **Apple Silicon CI** job.

Always: `verify_impaqt.sh` after every change; `bench_impaqt.sh` for anything touching
DBSCAN/paths or the assignment path.
