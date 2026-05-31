# Migration scripts

Ready-to-run helpers to finish wiring the two model repos to `dssatutils`.
Run them **in order**, reviewing the diff after each. All are idempotent and
back up every file they touch into a timestamped `.migration_backup_*` dir
inside the target repo.

| Script | What it does | Needs remote? |
|---|---|---|
| `00_create_remote.sh` | Create the **private** GitHub repo and push code + `v0.1.0` tag (gh or manual) | creates it |
| `01_local_install_check.sh` | Smoke-test that the package imports in Python and all R files source cleanly | no |
| `02_phase2_ml_phenology.sh` | Rewire **DSSAT ML Phenology Prediction**: remove duplicated scripts, add guarded `library(dssatutils)` | no¹ |
| `03_phase3_gridded.sh` | Rewire **DSSAT Gridded Run Tutorial**: swap `source()`/imports for the package, delete moved files, keep landcover, add pins | no¹ |

¹ Phases 2 & 3 default to installing from `alwinhopf/dssatutils@v0.1.0` (needs the
remote). To run them **before** the remote exists, install from the local path:

```bash
USE_LOCAL=1 ./migrate/02_phase2_ml_phenology.sh
USE_LOCAL=1 ./migrate/03_phase3_gridded.sh
```

## Recommended order
```bash
./migrate/01_local_install_check.sh          # verify package is sound
./migrate/00_create_remote.sh                # create private repo + push + tag
./migrate/02_phase2_ml_phenology.sh          # then review: git -C "../DSSAT ML Phenology Prediction" diff
./migrate/03_phase3_gridded.sh               # then review: git -C "../DSSAT Gridded Run Tutorial" diff
```

After Phase 2, in R: `renv::init(); renv::snapshot()` (that repo had no lockfile).

## Rollback
Each script copies originals into `<repo>/.migration_backup_<timestamp>/`.
To undo, restore from there or `git checkout -- <file>` (edits are staged, not committed).
