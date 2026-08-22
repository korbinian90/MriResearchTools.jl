# MriResearchTools.jl - issue list

Working notes from the architecture review of 2026-08-19..21. Uncommitted on
purpose. Cross-repository items and release ordering are in `issues-stack.md`
(X0..X7). Full evidence:
https://claude.ai/code/artifact/1dcb7a4a-4523-46fc-af23-7c5b0ecc3688

Measured on this machine unless marked *unverified*. State refers to branch
`claude/julia-repos-architecture-review-tj1zzz`.

This is the widest surface in the stack and the package every other one leans
on, which is why the unguarded parts of it are the biggest single risk in the
Julia half.

---

## Done on this branch

- **F1 misspelled keyword silently disabled global phase correction.**
  `ext/QSM_common.jl:9` called `romeo(...; correct_global=true)`; ROMEO reads
  only `:correctglobal`. Because ROMEO probes an untyped `keyargs...` with
  `haskey`, the wrong name was accepted and dropped, so global n*2pi correction
  had never run on `qsm_romeo_B0`, which is the path CLEARSWI's `qsm_contrast`
  uses. Verified empirically: misspelled output was identical to passing
  nothing. Honest caveat: on the bundled test data the correction is itself a
  no-op (the unwrapped median rounds to zero), so the fix restores intent with
  no visible output change here. Whether it changes results is data-dependent.
- **F2 partly fixed.** `ext/QSM_common.jl` was `include`d into both `QSMExt` and
  `QuantitativeSusceptibilityMappingTGVExt`; Julia itself reports
  `Method definition qsm_average(...) overwritten ... (check for duplicate calls to include)`.
  The backend-independent helpers (`qsm_mask_filled`, `weighted_average`,
  `laplacian_combine`) moved to `src/qsm_common.jl` and are defined once;
  `qsm_mask_filled` now also works with no backend loaded.
- **F9** `mask_from_voxelquality` and `romeovoxelquality` are `const`.
- **F16** `robustrescale!` no longer mutates a caller-supplied mask.
- **Five assertion-free test files now assert** (see below).
- **`src/citations.jl`** registers `:aspire` (carrying the patent notice),
  `:homogeneity`, `:laplacian`, `:rts`, `:phase_based_masking`, `:qsmxt`, plus
  `describe_input(path)` which adds NIfTI dimensions to the provenance record.
  Compat `ROMEO = "1.5"`; version 3.6.0.

## Open

### M1. `qsm_B0` and friends still collide on load order (F2 remainder)
The helpers were deduplicated, but `qsm_B0`, `qsm_romeo_B0`,
`qsm_laplacian_combine` and `qsm_average` are still defined with identical
signatures in both extensions, so whichever loads second wins. The package's own
`test/Project.toml` lists both backends. Papering over it is worse than naming
it: what is missing is a way for the caller to *select* a backend. Options are
dispatch on a backend type/singleton, or a keyword. That is an API decision, so
it is documented in the extension header and README rather than guessed at.
`ext/QSMExt.jl:7`, `ext/QuantitativeSusceptibilityMappingTGVExt.jl:7`.

### M2. F18 - `mcpc3ds_meepi` emits NaN where `mcpc3ds` does not
On the bundled test data: NaN in 13,938 of 639,846 voxels (2.2%), none from the
ordinary path on the same input. Not confined to background air; magnitude at
NaN voxels reaches 58% of the image maximum, so they sit inside the object.
Recorded as `@test_broken all(isfinite, corrected_me)` so a future fix trips the
test rather than passing unnoticed. Not fixed here because the cause is not
obvious from comparing the two code paths and guessing would be worse than
naming it. `src/mcpc3ds.jl:78`.

### M3. Abstract and untyped public API (F10)
Every CLI carries a `Dict{String,Any}` of settings; the option surface has no
type that can reject a wrong name or a wrong type. This is the structural cause
of F1: mistakes surface as silently different images rather than errors. Also
the thing that makes dispatch statically resolvable for X6.

### M4. Function barrier around NIfTI reading
Only worth doing when X6 is picked up, but it is the single biggest static-
compilation blocker on this side: `niread`/`readmag` accounted for 32 of the 68
verifier errors, and `byteswap` for 8 more, all from the element type and
endianness being decided by the file header at run time. Standard fix is to read
the header, then dispatch once into a concretely typed kernel. About half the
remaining work is inside NIfTI.jl, which is not ours: that is the honest limit
without upstream changes or a narrow typed reader of our own.

### M5. Base.pi shadowing (F16 remainder)
`src/niftihandling.jl:40`, `fix_ge_phase!` shadows `Base.pi` with a data-range
value, so the literal `2pi` in that scope means "full data range". Correct as
written, and one rename away from a subtle disaster.

### M6. Aqua, JET, property tests
None of the three. For the package everything else depends on, this is the gap
worth closing first (F14).

### M7. Unmeasured: type stability of the hot paths
`gaussiansmooth3d` (box filter, uses `DataStructures.CircularBuffer` in the hot
loop) and `robustmask` were never profiled for inference here. Worth an
`@code_warntype` / JET pass before assuming anything. *unverified*

---

## Test suite: what changed and what it now guarantees

The suite reported 80 passing assertions across 11 files. **Five files contained
zero assertions** and ran for 14 of the suite's 19 minutes: they called
functions and checked nothing, so they could only fail by throwing. That is why
F1 survived, the typo sat inside `qsm_romeo_B0`, exercised by a 7.5-minute test
that asserted nothing.

| File | Assertions | Invariant that carries the weight |
|---|---|---|
| `VSMbasedunwarping.jl` | 0 -> 9 | A zero voxel-shift map must be the identity, in either readout direction |
| `intensitycorrection.jl` | 0 -> 11 | Correction must make in-object signal more uniform: coefficient of variation has to drop. Plus `makehomogeneous` is exactly "divide by `getsensitivity`" |
| `mcpc3ds.jl` | 0 -> 18 | `combinewithPO`'s complex and PhaseMag methods compute the same sum, so they must agree |
| `qsm.jl` | 0 -> 9 | Finite, mask-confined, within a few ppm; estimators must correlate in-mask |
| `qsm_tgv.jl` | 0 -> 9 | Same, deliberately backend-agnostic because M1 means load order decides who answers |

80 -> 165 assertions, all executed.

**Two things deliberately not asserted**, because asserting them would have
meant asserting something untrue:
- `mcpc3ds_meepi`'s NaN output, recorded as `@test_broken` (M2) rather than
  accommodated.
- The correlation between supplied-mask and automask QSM runs. Those use
  genuinely different masks and QSM is strongly mask-sensitive: r = 0.37 on TGV,
  0.15 on QSM.jl. Any bound tight enough to mean something would really be
  asserting which backend loaded. The numbers are documented instead.

`Statistics` had to be added to `test/Project.toml`; ROMEO is reached as
`MriResearchTools.ROMEO` rather than added as a test dep.
