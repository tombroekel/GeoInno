# AI assistance in the GeoInno 2.0 Structural Diversity revision

## Methodological and code-review assistance

OpenAI ChatGPT, model GPT-5.6 Sol, assisted in October 2026 with review of the historical code, comparison with the intended NDS procedure, identification and discussion of sampling, RNG, undefined-value, and transformation issues, robustness-analysis design, and documentation planning. This assistance did not independently validate any scientific conclusion.

## Implementation assistance

OpenAI Codex assisted with implementation and tests on 26 October 2026. Model/version: not exposed by the Codex environment. The work was prepared on branch `structural-diversity-2.0` from base commit `e4819aeded4bb2112e987dd443a5f01030252043`. The published implementation and test commit is `57c8d4d90b567fb526c3dd62fe87c5c3b4cbd2e9`.

Codex modified `R/complexity_calculation.R` to retain the historical public estimator under `structural_diversity_legacy()` and add fixed-unique-node sampling, missing-value treatment, diagnostics, and reproducible RNG to the recommended `structural_diversity()`. It also added baseline and corrected-method tests, updated generated R help, README, NEWS, and package version metadata, and prepared comparison/audit material.

All AI-generated or AI-modified code is subject to human review and package tests. The package's human authorship and scientific responsibility are unchanged.
