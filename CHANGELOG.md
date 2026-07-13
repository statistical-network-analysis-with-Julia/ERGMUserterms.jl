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

- **Third-party terms are now first-class** (issue #1). ERGM.jl's term traits
  are a public, documented protocol (`ERGM.required_vertex_attributes`,
  `required_edge_attributes`, `requires_directed`, `requires_undirected`,
  `is_dyad_dependent`, `Networks.supports_missing`), and ERGM's formula
  validation and materialization act on the *declarations* rather than on its
  own term types. A custom term that declares its attributes is validated at
  `ERGMModel` construction exactly like a built-in one, instead of silently
  passing through and fitting an all-zero design column.
- **`validate_traits(term, net)`** — new, and run by `validate_term` unless
  `traits=false`. It exercises every declaration: that the declared
  vertex/edge attributes are the ones the term actually reads (attributes are
  perturbed and the statistic watched; reading an *undeclared vertex*
  attribute fails validation); that a `requires_directed`/`requires_undirected`
  declaration is enforced by `ERGMModel` on the direction twin of the network;
  that an `is_dyad_dependent = false` claim survives toggling every other dyad;
  that a `supports_missing = true` claim survives flipping the face value of a
  masked dyad; and that `ERGMModel` construction accepts the term.
- **`examples/MyTermPackage/`** — a package template for shipping third-party
  terms: Project.toml, module, and test suite for `ReciprocatedHomophily`, a
  term that declares all four traits, passes `validate_term`, and fits inside
  an `ERGMModel`. It is `include`d by this package's test suite, so it stays
  correct.
- `InteractionTerm` declares its two vertex attributes
  (`ERGM.required_vertex_attributes`), so a model naming attributes the
  network lacks now raises `ArgumentError` at construction. `WeightedEdges`
  deliberately declares no required edge attribute: it defaults the weight of
  edges that lack it, so it is well defined without the attribute.
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
