# Shared synthetic fixtures

`catalog.txt` has IDs 1–23. No observations were obtained from a real campaign.

- `split.cc`: two ten-node nearly complete graphs, each missing one edge. Clean values are `0.013*(i-j) + 0.0001*sin(i*j)` seconds, with S values scaled by 1.7 and quality 0.9. Contradictory ±10-second links join the groups and attach event 21; a zero-weight observation touches event 22; event 23 is absent. At `k_sigma=4` and minimum 30 edges, pruning leaves 44 edges in each group. Expected references are 3 and 13, chosen by surviving degree. Events 21–23 remain NaN.
- `split_tree.cc`: two triangles and a two-node tree connected by contradictory links. With minimum 3 edges, the triangles survive with references 1 and 4; the separated tree is rejected for insufficient redundancy.
- `median.cc`: the P observations from `split.cc` plus an S triangle. With exact size limit 2 and minimum 3 edges, P uses the finite graph bound while S uses exact uncertainty. This checks median weight scaling when the methods differ.

The Julia and native MATLAB regression suites read these same files. Probe studies generate additional graph structures at runtime with seeds 1–10.
