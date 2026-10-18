# MATLAB helper scripts

For synchronization, use `../run_ddsync.m` or the Julia implementation. `synchro_polished.m` is now a compatibility entry point into `ddsync.run`; it no longer contains a separate numerical implementation. Its legacy default input is `dt_long.cc`. Edit its `cfg` fields or migrate them to `ddsync_params.m`. The retired `STD_MIN_DIAGREL` clipping behavior is not available.

The following are dataset-specific preprocessing examples, not general catalog readers:

- `reindext_dtcc.m` reads `catalog_long.txt` and `dt_long.cc`, and writes **`catalog.txt` and `dt.cc`**. It assumes positive numeric IDs, retains the first duplicate catalog ID, and uses double-precision parsing. Do not use it for arbitrary large IDs, and do not overwrite the only copies of your inputs. Inspect the mapping and change filenames before running. Julia performs exact decimal-ID mapping automatically and writes `event_id_map.csv`.
- `reorder.m` moves column 7 of `catalog_long.txt` to the final column and writes `catalog_long_reordered.txt`. Use it only when column 7 is actually your event ID. It reformats numeric fields, so it is unsuitable for preserving arbitrary large integer IDs exactly.

Keep the catalog and both pair-header IDs consistent through any reindexing. Do not renumber the catalog alone after generating theta files.
