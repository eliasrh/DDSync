# DDSync (MATLAB + Julia)

This repository contains two parallel implementations of **DDSync**, a graph-based synchronization method for differential travel times. Both versions read a `dt.cc` file, solve per-(station,phase) graph synchronization problems, and write synchronized `dt_sync.cc` plus per-group `theta` outputs.

- **MATLAB implementation:** `DDSync/`
  - Mature, reference implementation.
  - Includes an example dataset (Spanish Springs; see below).
- **Julia implementation:** `DDSyncJulia/`
  - Generally faster runtime and lower overhead for large catalogs.
  - Offers a **dense vs. sparse storage option** to trade memory vs. speed for pass-2 output writing.

If you are new to DDSync, start with the README inside each subfolder:

- `DDSync/README.md` (MATLAB)
- `DDSyncJulia/README.md` (Julia)

## Which version should I use?

| Question | Suggested version |
| --- | --- |
| You already use MATLAB and want the original reference | MATLAB (`DDSync/`) |
| You need faster runtime and can use Julia 1.9+ | Julia (`DDSyncJulia/`) |
| You want example data included | MATLAB (`DDSync/`) |

Both versions aim to match **default parameter values** and **output formats**. If you find a mismatch, please report it so it can be fixed consistently.

## Example data (MATLAB only)

The MATLAB folder includes example `dt.cc` and `catalog.txt` from the Spanish Springs sequence (Trugman & Shearer, 2017: https://github.com/dttrugman/GrowClust). The event IDs have been made sequential and lightly filtered for demonstration. The Julia version does **not** include these files, but it reads the same formats.

## References

- Heimisson, E. R., Yu, Y. (2026, in preparation). *Graph-based denoising of differential travel-time observations with applications to pick reconstruction and path-difference tomography.*
- Trugman, D. T., & Shearer, P. M. (2017). *GrowClust: A hierarchical clustering algorithm for relative earthquake relocation, with application to the Spanish Springs and Sheldon, Nevada, earthquake sequences.* **Seismological Research Letters**, 88(2A), 379–391.


## License (non-commercial)

DDSync is licensed for **non-commercial** research and educational use only. Commercial or for-profit use (including internal commercial workflows) requires a separate written licensing agreement with the author (eliasrafn@hi.is). See `LICENSE` for full terms.
