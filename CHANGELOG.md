# Changelog

All notable changes to ERGMUserterms.jl are documented in this file. The
format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and the package adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] - Unreleased

Release driven by the 2026-07 expert-panel review: the custom-term contract
moves to ERGM.jl's state-independent add-direction `change_stat` convention,
and the validation harness now enforces it.

### Breaking

- **The `change_stat` contract changed to the add-direction convention.**
  Custom terms must return `g(y⁺ᵢⱼ) − g(y⁻ᵢⱼ)` — the statistic's change when
  edge (i, j) is present versus absent — *independent of the dyad's current
  state*. The old templates returned a toggle-signed value
  (`has_edge(net, i, j) ? -delta : delta`); every bundled example/template
  term has been rewritten, and the docs flag the old pattern as wrong.
  *Migration:* delete the `has_edge`-based sign flip from your
  `change_stat` methods and always return the add-direction value — the MH
  sampler negates it internally for removal proposals; re-run
  `change_stat_check`/`validate_term` afterwards.
- **The validation harness enforces state-independence.** `validate_term`,
  `change_stat_check`, and `consistency_check` brute-force
  `compute(present) − compute(absent)` and reject toggle-direction terms
  that previously passed. *Migration:* same as above — port your terms; a
  term that newly fails validation is following the 0.1 convention.
- **Minimum Julia raised to 1.12**; package UUID regenerated. *Migration:*
  upgrade Julia and re-resolve environments pinning the old UUID.

### Added

- Re-exports `compute`, `change_stat`, and `name` (imported from ERGM.jl),
  so term authors extend the ERGM generics directly with one `using`.
- `profile_term(term; sizes, density, n_iter)` is now implemented (the name
  was previously exported without a definition), profiling
  `compute`/`change_stat` cost across network sizes.
- `@ergm_term` now checks the generated definition (subtyping,
  `compute`/`change_stat`/`name` methods) and warns about missing pieces.
- Bundled terms declare `ERGM.is_dyad_dependent(...) = false`, opting out of
  the pseudo-likelihood caveat where appropriate; `WeightedEdges` gains a
  `default` weight keyword.

### Changed

- Bundled example terms (`WeightedEdges`, `DyadCovTerm`, `InteractionTerm`)
  read undirected edge attributes by canonical `(min, max)` key with
  `get(..., default)` fallbacks.
- Documentation (term interface, templates, validation) rewritten around
  the add-direction convention, including an explicit migration warning.

### Fixed

- `WeightedEdges.compute` sums over the network's actual edges with a
  default weight, so `compute` and `change_stat` stay consistent when the
  sampler adds edges without a weight attribute.
- The validation harness snapshots and restores dyad edge attributes around
  its toggle test, preventing spurious failures for attribute-reading
  terms.

## [0.1.0] - 2026-02-09

Initial release: `@ergm_term` scaffolding, example terms, and the term
validation/benchmark harness.
