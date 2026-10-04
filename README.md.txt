# Navigator's Log R&D: rocky-planet formation runs

Status: DRAFT [CL], 4 October 2026. Nothing here is frozen or certified yet.

Part of the Tunable Rocky Planet Water Model (Navigator's Log R&D, Christopher Blake Head, ORCID 0009-0004-2308-6051). This repository holds the N-body accretion runs behind decision D10: the model grows its own rocky planets instead of importing published outcomes. The specification and the Stage 2 working document (process inventory, decisions D10 to D12, the P-D12 preregistration) live in Google Docs and are the source of record for every decision.

Decision markers: [CBH] marks a decision by Christopher Blake Head; [CL] marks a decision or proposal that originated with Claude and stands until CBH confirms it.

## Pinned toolchain
- REBOUND 4.4.10 (C version), integrator TRACE.
- Fragmentation module by Childs and Steffen (2022, MNRAS 511:1848), github.com/annacrnn/rebound_fragmentation, commit d67ab2b, used unmodified.
- scripts/setup.sh fetches both, builds src/formation_run.c, and writes SHA-256 hashes to build/build_hashes.txt.

Known issue [CL]: in the module's hit-and-run branch, the interacting fraction beta is computed without density, so it is not the pure fraction its paper defines. P-D12 uses the module unmodified; any fix would be a new instrument and needs a new preregistration.

## What runs
- src/formation_run.c: 154-body Chambers (2013) bimodal disk (14 embryos of 0.093 and 140 planetesimals of 0.0093 Earth masses, 0.3 to 2.0 AU), Jupiter and Saturn on current orbits, 6-day step, density 3 g/cm3, minimum fragment mass half a planetesimal. Arguments: seed, expansion factor, simulated-time cap, wall-clock limit. It saves and resumes on its own, so long runs are split into sessions.
- Two inputs are inferred, not read from the source, and are flagged in the code: the body masses (from the published total) and the surface density starting from zero at 0.3 AU.

## Runs planned
1. B0 (exploratory): one real-size and one 3x run from seed 900001, one 6-hour session each, via the "B0 benchmark" workflow. It measures cost and sets the P-D12 matched-stage fraction. It is never evidence about inflation, and its seed is never reused.
2. P-D12 (confirmatory, preregistered): 20 paired runs (real size and 3x, same seed per pair), 10 extra real-size control runs, and one same-seed repeat. The analysis script is frozen and hashed, and the preregistration must pass its certification check, before any P-D12 run starts.

## Publish either way
The P-D12 verdict, accepted, not accepted, or inconclusive, will be reported here at equal prominence.

## License
To be decided by CBH. Because formation_run.c compiles the GPL-3.0 fragmentation module into the program, GPL-3.0 is the straightforward choice for this repository [CL].
