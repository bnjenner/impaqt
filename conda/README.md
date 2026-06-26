# Conda recipe (bioconda)

Draft recipe to distribute impaqt via [bioconda](https://bioconda.github.io), so
users can `conda install -c bioconda impaqt` and get a prebuilt binary for
linux-64, linux-aarch64, osx-64, and osx-arm64 (bioconda's CI builds each).

## Files
- `meta.yaml` — package metadata, dependencies, source tarball + checksum, test.
- `build.sh` — build steps; uses the **conda-provided bamtools** and skips the unit
  tests (`-DUSE_SYSTEM_BAMTOOLS=ON -DIMPAQT_BUILD_TESTS=OFF`) so nothing is fetched
  from the network during the build.

## Before submitting
1. **Release**: pinned to `v1.2.0` with its tarball `sha256` already filled in
   `meta.yaml`. To re-point at a future release, bump `version` and recompute:
   ```
   curl -sL https://github.com/bnjenner/impaqt/archive/refs/tags/v<VER>.tar.gz | sha256sum
   ```
2. (Recommended) test the recipe locally:
   ```
   conda install -n base conda-build
   conda build conda/
   ```

## Submitting
bioconda recipes live in their own repo, not here:
1. Fork https://github.com/bioconda/bioconda-recipes
2. Copy these files to `recipes/impaqt/{meta.yaml,build.sh}`
3. Open a PR; bioconda CI builds + tests all platforms. On merge, it's published to
   the bioconda channel.

This `conda/` folder just keeps the recipe under version control alongside the code.
