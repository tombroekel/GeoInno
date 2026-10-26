# GeoInno 2.0.0

This is a breaking methodological release for Structural Diversity. The historical implementation remains available through `structural_diversity_legacy()` specifically to reproduce analyses based on earlier GeoInno versions.

- `structural_diversity()` now samples up to 125 starts and walks until 150 distinct vertices have been visited, replacing the historical fixed 200-step walk. Components no larger than the requested node count are used in full.
- Random-number handling now uses a top-level `seed` and reproducible streams across sequential and parallel workers. The legacy estimator retains its historical seed resets.
- Undefined sample-level iNDS values are missing rather than substituted with zero. Valid sample-level values are averaged; diagnostics report how many samples were valid.
- Raw NDS is no longer capped at 1. Both `Structural_raw` and its `Structural_neglog` transformation are returned; negative transformed values can occur when raw NDS exceeds 1.
- The historical `NDS.intern()` and `nds()` helpers remain numerically unchanged and are documented as legacy helpers.
- The methodological review and implementation were AI-assisted; see `AI_ASSISTANCE.md` for the scope and model information. Human authorship and responsibility remain unchanged.
