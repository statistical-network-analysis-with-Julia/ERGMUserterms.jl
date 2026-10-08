# Changelog

All notable changes to ERGMUserterms.jl are documented in this file. The
format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and the package adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] - Unreleased

First public release of ERGMUserterms.jl, the Julia counterpart of R's
`ergm.userterms`: an interface, a validation harness, bundled example terms
and a copyable package template for writing custom ERGM.jl terms. The
bundled terms and the harness's brute-force reference are checked against R
ergm 4.12 by a provenanced fixture. The changes below are relative to
0.1.0, which was not publicly released.

**Dependency renamed:** the foundation package is now `NetworkCore` (developed as `Networks`); write `using NetworkCore` where code said `using Networks`. Types and functions keep their names.

### Highlights

- **Add-direction change statistics.** `change_stat(term, net, i, j)` is
  `g(y⁺ᵢⱼ) − g(y⁻ᵢⱼ)`, independent of the dyad's current state — ERGM.jl's
  convention — and the harness enforces it.
- **A harness that checks every dyad.** `validate_term` compares the change
  statistic with a brute-force recomputation on every dyad of the network,
  checks the statistic of the empty network, and tests each trait
  declaration, including dyad-independence and the survival of edge
  attributes under MCMC.
- **Custom terms are first-class in ERGM.jl.** A term declares what it
  needs through ERGM.jl's public term traits and is validated at model
  construction exactly like a built-in term.
- **`WeightedEdges` is `edgecov(W)` under MCMC**, not only under MPLE.
- **Reproducible validation.** Every random draw comes from an `rng`
  keyword, and every failure message carries the literal that replays it.

### Breaking

- **The `change_stat` contract is the add-direction convention.** The 0.1
  templates returned a toggle-signed value (`has_edge(net, i, j) ? -Δ : Δ`);
  every bundled term is rewritten, and `validate_term`, `change_stat_check`
  and `consistency_check` reject toggle-direction terms. *Migration:* delete
  the `has_edge` sign flip and always return the add-direction value (the
  sampler negates it for removals). The same applies when porting a statnet
  C changestat: drop the `edgestate` sign.
- **Every dyad is checked by default.** `validate_term`,
  `change_stat_check` and `test_term` take `n_tests=nothing` (was 10 / 10 /
  100 uniformly random dyads): every dyad of a network with up to 5000
  dyads, a stratified sample of 5000 beyond. An integer `n_tests` is now a
  stratified sample of that many *distinct* dyads (strata: the dyad's
  state, the reverse arc's state, shared partners). `consistency_check`
  defaults to `exhaustive=true`. A mutual-type term that was wrong on half
  of the reciprocated dyads passed the old default on 136 of 200 sparse
  20-vertex networks; it now fails on all 200. *Migration:* pass
  `n_tests=k` / `exhaustive=false` where the exhaustive run is too slow.
- **`validate_term` checks `g(∅)`.** `compute` on the empty network plus
  the change statistics summed while the edges are added one at a time must
  equal `compute` on the network. A `compute` that special-cases the empty
  network to 0 for a statistic that is not 0 there (an isolates count) now
  fails; the guide no longer recommends that shortcut.
- **A mis-declared `is_dyad_dependent = false` always fails
  `validate_traits`.** The change statistics are compared on the empty, the
  observed and the complete network, and every other dyad is toggled for
  every source dyad (all pairs up to about 500 dyads; the reverse arc and
  incident dyads first beyond that). The old check toggled around four
  random dyads and missed a mutual-type term on 18 of 100 seeds. ERGM.jl now
  refuses such a term when the model is built; `validate_traits` reports the
  mismatch once, not again as an `ERGMModel` construction failure.
- **A term that reads an edge attribute live fails `validate_traits`** on a
  network carrying edge attributes: ERGM.jl's samplers delete edge
  attributes when they toggle a dyad. *Migration:* hold the attribute as a
  matrix (see `WeightedEdges` below).
- **`DyadCovTerm` refuses a matrix of the wrong size** (and a non-square
  one) with an `ArgumentError` from the constructor, `compute`,
  `change_stat` and `ERGMModel` construction. Entries outside the matrix
  used to count 0, and a 4×4 matrix on 16 vertices fitted to `-Inf`.
- **Networks ERGM.jl cannot fit are refused up front** by every entry point
  that takes a network, with an `ArgumentError`: fewer than 2 vertices, a
  two-mode network, or a network containing a self-loop (ERGM.jl's own
  messages for the last two).
- **No verdict without evidence.** An integer `n_tests` and `n_iter` must
  be at least 1; `test_term` refuses `n_vertices < 2` and a `density`
  outside `[0, 1]`.
- **A declared attribute that is constant on the validation network fails
  `validate_traits`**: whether the term reads it cannot be decided there. A
  ✓ line is printed only for a check that ran.
- Minimum Julia is 1.12; the package UUID was regenerated.

### Added

- **`validate_traits(term, net)`**, run by `validate_term` unless
  `traits=false`: the declared attributes are the ones the term reads, and
  are complete; direction requirements are enforced by `ERGMModel`; the
  dyad-independence and `supports_missing` claims hold; edge attributes
  survive the sampler; `ERGMModel` accepts the term.
- **Term traits.** Terms declare `ERGM.required_vertex_attributes`,
  `required_edge_attributes`, `requires_directed`, `requires_undirected`,
  `is_dyad_dependent` and `NetworkCore.supports_missing`; the bundled terms
  declare theirs.
- **`WeightedEdges` snapshots its weights at model construction**
  (`ERGM.Extension.materialize`): inside a model it is `edgecov(W)` in MPLE, MCMLE,
  simulation and goodness of fit. It gains a `default` weight keyword and
  declares no required edge attribute.
- **`examples/MyTermPackage/`**, a package template for third-party terms
  (`ReciprocatedHomophily`), tested by this package's suite and in CI.
- **`rng::AbstractRNG` on every validator and benchmark**, with a replay
  literal (`reproduce with rng=Xoshiro(0x…)`) in every failure message.
- `test_term` and `profile_term` take `directed=` and `vertex_attributes=`,
  so attribute-declaring and undirected-only terms can be tested and
  profiled; `test_term` honours `n_tests`.
- `profile_term` (exported but undefined in 0.1) is implemented.
- `@ergm_term` checks the definition and warns about missing methods.
- A term whose methods were defined after a bare `using` (local functions
  ERGM.jl never calls) gets a message naming the missing
  `import ERGM: name, compute, change_stat`.
- `compute`, `change_stat` and `name` are re-exported; they are ERGM.jl's
  generics.
- Masked dyads: a term without `supports_missing` is validated at their
  face value, and the harness says so.
- **Validation against R**: `test/fixtures/userterms_examples.toml` (ergm
  4.12.0) pins every bundled term's statistic and change statistics on
  every dyad, and the harness's brute-force reference, against the statnet
  terms they equal (`nodecov`, `edges`, `edgecov`, `mutual(same=)`), on a
  directed network and its undirected projection; and the fitted
  `edgecov(W)` model, which `WeightedEdges` reproduces by MPLE and MCMLE.
- Allocation and inference pins: every bundled term, the template's term
  and the `WeightedEdges` snapshot are 0 B per `change_stat`/`compute` and
  infer `Float64` (`benchmark/regression_tests.jl`, run in CI); a
  BenchmarkTools suite checks the per-dyad cost does not grow with n.
- Every export has a docstring with a runnable example, and every README
  and guide example is executed by the test suite, which also checks the
  outputs the pages quote.

### Changed

- The harness sizes its random dyad budget with ERGM.jl's
  `ERGM.Extension.n_observed_dyads` (the observed dyads, masked ones
  excluded) instead of a private copy, and refuses a two-mode or looped
  network with `ERGM.Extension.require_supported_network`, the check
  `ERGMModel` runs.
- A term that snapshots an attribute at model construction adds a method to
  `ERGM.Extension.materialize`, ERGM.jl's extension API, where development
  versions extended the then-`public` `ERGM._materialize` (the templates
  guide and `WeightedEdges`/`DyadCovTerm` show it).
- The workflow testset skips the CI layout step, with a message, when no
  sibling checkout is beside the package (a lone checkout or a registry
  install) instead of failing.

- Bundled terms read attributes per element and infer `Float64`
  (`WeightedEdges` and `InteractionTerm` inferred `Any` and allocated on
  every call).
- `consistency_check` visits each unordered dyad once on an undirected
  network; its random mode draws distinct off-diagonal dyads.
- The documentation is rewritten around the add-direction convention, the
  term traits, edge-attribute snapshots and one import idiom
  (`import ERGM: name, compute, change_stat`); installation follows the
  ecosystem workspace recipe.
- A term that reads an edge attribute live is now refused by ERGM.jl's
  samplers (MCMLE, `simulate_ergm`, `gof`, the MPLE bootstrap) instead of
  being simulated under its all-default model; the documentation says so.
- ERGM.jl's `fit_ergm` defaults to `method=:auto` (R's rule): a formula
  with a user term that declares nothing is fitted by MCMLE, since an
  undeclared term counts as dyad-dependent. Declare
  `ERGM.is_dyad_dependent(::MyTerm) = false` for a covariate-only term
  (the MPLE is then exact), or pass `method=:mple`.
- The harness works with ERGM.jl 0.2's fallbacks, which are
  `ArgumentError`s typed on the term alone: a `change_stat(::MyTerm, net,
  i, j)` with unannotated `i, j` now fits end to end.

### Fixed

- `WeightedEdges` gave wrong MCMLE, `simulate_ergm` and `gof` results: the
  sampler's toggles deleted the stored weights, so the chain targeted the
  all-default model (MCMLE −1.266 with `converged = true` against an exact
  −1.236 on a 20-vertex example).
- `WeightedEdges.compute` sums over the network's edges with the default
  for unweighted ones, consistently with `change_stat`.
- The harness restores a dyad's edge attributes around its toggles, so
  attribute-reading terms are not reported inconsistent.
- Documentation examples that could not run as written, and output
  comments the code did not produce.

### Known limitations

- R's C template (`changestats.users.c`) and its `mindegree` example term
  are not ported; `examples/MyTermPackage/` is the Julia counterpart of the
  template.
- A pass is a statement about the network validated on: a change statistic
  wrong only in a configuration that network lacks is not seen.
- Above 5000 dyads the per-dyad checks use a stratified sample, and the
  dyad-independence check toggles all pairs of dyads only up to about 500
  dyads (22 vertices directed, 32 undirected).
- ERGM.jl's samplers do not preserve edge attributes. `validate_traits`
  fails a term that reads one live, and ERGM.jl's samplers refuse such a
  term with an `ArgumentError` (the MPLE accepts it). ERGM.jl trusts a
  term's `is_dyad_dependent` declaration beyond its probe at model
  construction.
- The harness checks computation, not identifiability, non-degeneracy or
  model fit.

## [0.1.0] - 2026-02-09

Initial development version (not publicly released): `@ergm_term` scaffolding, example terms, and the term
validation/benchmark harness.
