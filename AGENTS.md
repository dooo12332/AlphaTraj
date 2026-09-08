# AGENTS.md

AlphaTraj ("pocket analysis"): research CLI/Python tool that analyzes dynamics of protein sub-pockets from MD trajectories. No tests, no CI, no linter, no packaging, no requirements file. User guide: `tutorial/AlphaTrajTutorial_ENG.md` (authoritative for parameter semantics); `README.md`.

## Layout & architecture

- `pocket_analysis.py` — the entire tool (single 1296-line module, no package). Top-level `argparse` setup runs on import, but execution only under `if __name__ == '__main__'` (line ~1215), so it can be `import`ed as a library.
- `alphaspace2/` — a **modified copy of AlphaSpace2.0** (do NOT swap in the upstream RedesignScience package; it is incompatible). Runtime deps: `numpy`, `scipy`, `mdtraj`, `tqdm`. `alphaspace2/data/*.dat` are required and must stay alongside the package.
- Main classes (all in `pocket_analysis.py`): `PocketsAnalysis` (:477, main entry), `PocketsInfo` (:77, aggregated per-pocket results), `SnapShots` (:17, per-frame), `ModelGroup` (:672, multi-model pocket comparison), `PAHelper`/`MGHelper` (writers), factory `GetPA` (:1037), `WriteFiles` (:1101), `ParserMask` (:1024).
- Unit quirk: `mdtraj` loads nm; internal geometry math multiplies by 10 to work in Å (e.g. `_CalcBoxEdge` :42). SASA calls in `alphaspace2/functions.py:getSASA` branch on the `mdtraj._geometry._sasa` signature (1.10.0 vs 1.10.2+); do not "fix" the call to a single signature.

## Running

```sh
python pocket_analysis.py            # prints help (auto `-h` when no args)
python pocket_analysis.py --top T --traj A.nc,B.nc --box ...   # CLI mode (single model)
python pocket_analysis.py --config config.ini                  # ini mode (single or multi)
python pocket_analysis.py -v,--version
```

- Must run from the repo root (or add it to `sys.path`): the code does `import alphaspace2`.
- There is no example trajectory in the repo; an analysis run needs a user-supplied, pre-stripped topology + trajectory (solvent/ions removed, protein imaged/centered — see tutorial §Trajectory Preprocessing). Expect long runtimes; the per-frame AlphaSpace pass is the bottleneck. Prefer `--box`/frame extraction and save results with `pickle=<file>` so `unpickle=<file>` later skips `Analysis`.
- No way to run "a test" — there are none.

## Two input modes

- **CLI flags** (single model only). Masks are space-separated on the CLI.
- **`--config` ini** via `configparser`: sections `[GENERAL]` + `[MODEL0]`, `[MODEL1]`, … (`MODEL<n>` loop caps at 100). `[GENERAL]` needs `mode=single|multi`. `multi` compares several systems: requires `pock_pairs=<outfile>` in `[GENERAL]` and an `align_mask` in every model section (same residue ranges per model; aligned on backbone CA/C/N only). This mode is only reachable via ini — CLI cannot express it.
- Library use mirrors ini keys as a `dict[str,str]`: `pa = PA.GetPA(param_dict, PA.ParserGeneral({}))`, then `pa.Analysis(start, stop, offset)`; reuse with `{'unpickle': './raw.dat'}`. See tutorial end for a full example.

## Mask/box semantics (easy to get wrong)

- `rec_mask`/`lig_mask`/`align_mask`: residue ranges in **pairs**, 1-based, start **inclusive**, stop **exclusive**, comma-separated. `30,51,60,61,80,102` = residues 30–50, 60, 80–101. Always an even count (a single residue needs two numbers). If neither mask given, everything is receptor; if only `lig_mask` given, receptor = all atoms minus ligand.
- `--box atom1[,atom2…] lenX,lenY,lenZ` — one or more boxes; atom ids are 1-based atom serials. Box centers follow the chosen atoms' geometric center per frame.
- When `lig_mask` is set, occupancy/contact analysis is enabled; without it occupancy outputs are skipped.
- Frame args: `frame_start` inclusive, `frame_stop` exclusive, `frame_offset` stride.

## Gotchas

- **CLI boolean flags are broken for disabling**: argparse declares `type=bool`, so `--out_xxx false` evaluates truthy. To turn outputs off, use the `--config` ini (values go through `str2bool`) or the Python API. This also means tutorial mentions of `--out_corr` don't exist in argparse (the real flag is `--out_coex`); tutorial text is not always in sync with the code.
- Only `pock_info/summary.dat` references the internal pocket **id**; every other output file references the **rank**. Pocket ids in files are rank unless it's summary.dat.
- Outputs land under `--out/pock_info/`: `summary.dat`, `best_conf/` (`out.dat`, `pdb/`, `vol_*_v_time.dat`), plus optional `pock_data/` (`.npy`/`.txt`), `main_group/`, coexistence/correlation/group/MTP matrix `.dat` files.
- Score cutoffs: score ≥ 20 = "main" sub-pocket, 10–20 = "minor", < 10 = invalid. See `PocketsInfo._CalcPocketScore` / `WriteSummaryInfo` for the conventions.
- Existing commits are descriptive, in Chinese for non-`mdtraj` fixes; match repo style when committing.
