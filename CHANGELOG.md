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
  `change_stat_check`/`validate_term` afterwards. If you are porting a
  statnet changestat, note that the toggle-direction idiom is R ergm's own C
  convention (`C_CHANGESTAT_FN(…, Rboolean edgestate)`, `CHANGE_STAT[0] +=
  edgestate ? -1 : 1`): drop the `edgestate` sign and return the
  add-direction value.
- **Two-mode networks are refused.** `validate_term`, `validate_traits`,
  `change_stat_check`, `consistency_check` and `benchmark_term` throw an
  `ArgumentError` (carrying ERGM.jl's one-mode-only message) on a
  bipartite-flagged network (`network(n; bipartite=k)`). Previously they
  ran: their dyad draws range over every off-diagonal pair, and on a
  two-mode network the within-mode pairs cannot hold an edge (`add_edge!`
  returns `false`), so the brute-force reference read 0 there and a
  **correct** term was reported inconsistent (`Inconsistency at dyad (1,3):
  change_stat=4.0, brute-force=0.0`). ERGM.jl cannot fit such a network
  anyway. *Migration:* validate on a one-mode network.
- **The validation harness enforces state-independence.** `validate_term`,
  `change_stat_check`, and `consistency_check` brute-force
  `compute(present) − compute(absent)` and reject toggle-direction terms
  that previously passed. *Migration:* same as above — port your terms; a
  term that newly fails validation is following the 0.1 convention.
- **Minimum Julia raised to 1.12**; package UUID regenerated. *Migration:*
  upgrade Julia and re-resolve environments pinning the old UUID.
- **No verdict without evidence.** `validate_term`, `change_stat_check` and
  `test_term` require `n_tests >= 1` (`n_tests=0` used to mean "skip the
  only numeric check" and returned `true` — for a toggle-direction term
  too); `validate_term`, `validate_traits`, `change_stat_check`,
  `consistency_check` and `benchmark_term` throw an `ArgumentError` on a
  network with fewer than 2 vertices (the validators used to return `true`
  having checked no dyad; `benchmark_term` already refused); `test_term`
  refuses `n_vertices < 2` and a `density` outside `[0, 1]` (it printed
  `All tests PASSED` for `n_vertices=1`, `density=1.5` and `density=-0.2`,
  the change-statistic block silently skipped). A toggle-direction term can
  no longer pass any validator under any accepted keyword combination (a
  testset enumerates them). *Migration:* pass `n_tests >= 1` and a network
  with at least 2 vertices.
- **Networks containing a self-loop are refused up front**, like two-mode
  ones: every network-taking entry point throws an `ArgumentError` carrying
  ERGM.jl's own self-loop message (ERGM.jl models the off-diagonal dyads
  only; `rem_edge!(net, v, v)` or `loops=false` is the fix it names) — since
  the reconciliation round the message IS ERGM's, from its now-`public`
  `_refuse_self_loops` / `_refuse_two_mode`, prefixed with the entry point
  (`"validate_term validates terms on networks ERGMModel accepts only, and
  ERGMModel refuses this one: …"`) rather than a copy mirrored by hand. Before,
  such a network was validated dyad by dyad with `compute` counting the loop
  (every check green), then rejected at the `ERGMModel` step in words that
  attributed the network's defect to the term — and `traits=false`,
  `change_stat_check` and `consistency_check` returned `true` outright. A
  `loops=true` network with no loop present is unaffected.
- **A declared attribute that is constant on the validation network fails
  `validate_traits`** (and so `validate_term`): perturbing a constant
  attribute changes nothing, so whether the term reads it is undecidable
  there, and the network is not one "the term is meant for". Before, the
  check was skipped and the ✓ line `declared attributes match the ones the
  term reads` printed anyway — an *undeclared*-attribute term (the silent
  all-zero-column trap the check exists for) passed on such a network with a
  positive finding that was never made. An undeclared constant attribute now
  warns (`… is constant on this network, so whether the term reads it is
  undecidable; validate on a network where it varies`) and the ✓ is replaced
  by a line naming what could not be tested; the verdict stands.
  *Migration:* validate on a network where every declared attribute varies.

### Added

- `profile_term` takes `directed::Bool=true` and
  `vertex_attributes::AbstractDict{Symbol}=Dict{Symbol,Any}()`, the two
  keywords `test_term` gained, and builds its networks the same way. An
  undirected-only term whose `change_stat` refuses a directed network could
  not be profiled at all (it raised inside `benchmark_term`), and an
  attribute-declaring term was timed on its attribute-absent fallback branch
  rather than its real cost. A `ProbeTerm` testset asserts every generated
  network has the requested kind.
- **Output comments are gated.** The docs-execution testset now captures
  each page's stdout and `@info` records and requires every numeric
  `compute() returns: X` / `Works on empty|complete network: X` comment on
  the page to be a value some run on that page produced — the four numbers
  fixed below could not have survived it (a negative run with `30.0`
  reinstated fails with `validation.md:83 quotes \`compute() returns:
  30.0\` but no run on the page produced that value`).
- **The bare-`using` mistake gets an actionable message.** When no
  `compute`/`change_stat`/`name` method for the term reaches the shared
  generic — what a method defined after `using ERGM` without `import ERGM:
  …` looks like from inside the harness — `validate_term`'s failure lines
  and `@ergm_term`'s definition check append `no compute method for MyTerm
  reaches the shared generic; if you defined one after a bare \`using\`, it
  is a local function ERGM.jl never calls: add \`import ERGM: name, compute,
  change_stat\` before the definitions`. Before, the author read only
  ERGM's fallback `compute() not implemented for MyTerm` while looking at
  the `compute` method they had just written. A term without a `name`
  method is reported as using ERGM's fallback name (the run still passes).
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
- **`rng::AbstractRNG` keyword on every validator and benchmark** (panel
  2026-09, item 6): `validate_term`, `validate_traits`, `change_stat_check`,
  `consistency_check`, `test_term`, `benchmark_term` and `profile_term` all
  take `rng=Random.default_rng()` and draw every random dyad and every random
  test network from it — nothing from the global RNG. Two calls with equal
  rngs visit the same dyads and return the same verdict, at any thread count.
- **Every failure message reports how to replay it.** Each harness call
  snapshots its rng state on entry and every failure `@warn` ends with
  `reproduce with rng=Xoshiro(0x…, 0x…, 0x…, 0x…) passed to <entry point>`;
  passing that literal back visits the identical dyad sequence, failing dyad
  included. `test_term` prints the literal in its report header (its
  failures may sit in the network it generated). `consistency_check` gains
  `verbose=false` to report the failing dyad the same way.
- `n_tests` keyword on `validate_term` (default 10, capped at the number of
  dyads): the number of random dyads for the `change_stat` type check and
  for the forwarded `change_stat_check`, replacing the hard-wired 10 and 5.
- **`test_term` can test attribute-declaring and undirected-only terms:**
  new keywords `directed=true` and `vertex_attributes=Dict{Symbol,Any}()`
  (`attr => v -> value` or `attr => vector`), applied to the random, the
  empty and the complete network alike. Before, it could build only a bare
  directed network, so the package's own template term
  (`ReciprocatedHomophily`, which declares `:group`) and any
  `ERGM.requires_undirected` term "failed" the tutorial's step 3 with
  `term declares required vertex attribute :group, which the validation
  network does not have` / `validated on a directed network`. Those
  messages are still what you get without the keywords — the docstring,
  `validation.md`, `getting_started.md` and the README now say so.
- **Masked dyads are announced.** For a term with `supports_missing ==
  false` on a network with masked dyads, `validate_term`,
  `change_stat_check` and `consistency_check(…; verbose=true)` log one
  `[ Info: net has k masked dyads; the term declares supports_missing =
  false, so they are validated at their face value (ERGM.jl applies its
  missing= policy at estimation time)]` line — the missing-data contract's
  requirement that a face-value number never be handed out silently. The
  verdict is unchanged.
- **Inference pins on every bundled term.** `WeightedEdges` and
  `InteractionTerm` had `change_stat`/`compute` inferring `Any` (values
  pulled from the untyped `Dict{…,Any}` attribute stores with a bare
  `get`), which cost 96 B / 208 B of boxing per call and, downstream, made
  ERGM's statically typed change-statistic tuple infer `Tuple{Float64, Any}`
  for every model containing them (a dynamic convert on every MPLE row and
  MH step). They — and the template's `ReciprocatedHomophily` — now read per
  element with `get_vertex_attribute(net, attr, v)` /
  `get_edge_attribute(net, attr, i, j)` (which return `nothing` when absent,
  canonicalise the undirected key and allocate nothing) and assert
  `::Float64`. The allocation pins in `benchmark/regression_tests.jl` and the
  "Allocation regressions" testset are **0 B** for all six terms (`compute`
  too, on directed and undirected networks; an earlier 0.2.0 draft pinned
  the measured ≤ 96 / ≤ 208 / ≤ 80 B ceilings instead), and they gain
  `Base.return_types(change_stat, …) == [Float64]` /
  `Base.return_types(compute, …) == [Float64]` so a shipped template can
  never infer `Any` again. README "Benchmarking", `benchmarking.md` and
  CLAUDE.md "Performance pins" describe the new numbers.
- The package template's three `fit_ergm` calls (`test/runtests.jl`,
  `README.md`, the `ReciprocatedHomophily` docstring) pass `rng=Xoshiro(1)`,
  so the skeleton third-party authors copy models the rng contract end to
  end (MPLE is deterministic; `method=:mcmle` is not).
- `validate_traits` reports a declared vertex attribute that is set on only
  some vertices (ERGM.jl refuses partial attributes at model construction,
  as statnet refuses NA), naming the vertices without a value; its final
  `ERGMModel` step prints ERGM's own message rather than an
  `ArgumentError(...)` wrapper.
- `benchmark_term` refuses `n_iter < 1` and a network with fewer than two
  vertices with an `ArgumentError` instead of a `DivideError`/`BoundsError`.
- **Golden fixture against statnet** — `test/fixtures/userterms_examples.toml`,
  generated by the checked-in `test/fixtures/r/userterms_examples.R` (ergm
  4.12.0 / R 4.6.1, `[provenance]` block, loaded with `Networks.load_golden`),
  pins every bundled term's statistic and its add-direction change statistic
  on **every dyad** — and the validation harness's own brute-force reference
  (`_brute_change_stat`) — against the statnet terms they are identical to:
  `ExampleTerm` = `nodecov("id")`, `TemplateTerm(1.0)` = `edges`,
  `WeightedEdges()` = `edgecov(W)`, `DyadCovTerm(cov)` = `edgecov(cov)`,
  `InteractionTerm(:a,:b)` = `edgecov(M)`, and the template's
  `ReciprocatedHomophily(:group)` = `mutual(same="group", diff=FALSE)`, on an
  8-vertex directed network and on its undirected projection (which pins the
  `(min,max)` canonical-key branches). Tolerance `1e-9`: deterministic
  arithmetic on a fixed graph. R's coefficient labels are recorded but not
  claimed — user terms keep their own names.
- **`benchmark/` suite and allocation pins.** `benchmark/benchmarks.jl`
  (BenchmarkTools `SUITE`: `compute` and `change_stat` of every bundled term
  at n = 500 and 2000 with constant mean degree, `BENCHJL` rows for the site's
  `tools/run_benchmarks.jl`, and a ≤ 3× per-dyad scaling assertion) and
  `benchmark/regression_tests.jl`, mirrored by the testset "Allocation
  regressions": every bundled term and the template's `ReciprocatedHomophily`
  is pinned at **0 B** per `change_stat` and per `compute` (on directed and
  undirected networks) and at `Base.return_types(…) == [Float64]` for both
  methods (see "Inference pins" below for why). `_other_dyads`' 250-dyad
  cap is pinned on a 30-vertex network.
- **A runnable example on every export, executed by the test suite.** Every
  exported name — the seven validators/benchmarks, the five bundled terms,
  `AbstractUserTerm`, `@ergm_term` (whose example is now a real term,
  `ScaledEdges`, with a numeric `compute`/`change_stat` that `validate_term`
  accepts), `term_signature`/`term_documentation`, and the docstrings this
  package attaches to the re-exported generics `name`/`compute`/`change_stat`
  — carries a self-contained ```julia block that does its own `using` and
  ends in an assertion-style comment. The testset "Every exported docstring
  carries a runnable example" walks `names(ERGMUserterms)` and runs every
  block in a fresh module (19 blocks; `>= 18` pinned).
- **README/docs execution gate in the test suite.** The testset "README and
  docs pages execute without warnings" is the CI-runnable form of the site's
  `tools/check_snippets.jl`: it extracts the ```julia blocks of `README.md`
  and every `docs/src/**/*.md` with the checker's skip rules, runs each
  file's blocks in order in one fresh module (REPL soft scope, temp dir,
  stdout swallowed), and fails on any exception **and on any `@warn` whose
  `_module === ERGMUserterms`** — so an example that "runs" while
  `validate_term` returns `false` fails the suite instead of printing a
  warning. Every page claim that used to be a bare `validate_term(…)` call is
  followed by `@assert valid`. It also greps every doc line for
  `using .*\bNetwork\b` and for bare `rand(1:…)` draws, and requires each
  page that validates to pass `rng=Xoshiro(…)`.
- README and `docs/src/index.md` gain a **"Not implemented (vs R
  ergm.userterms)"** section: R's C `changestats.users` template and its
  `mindegree` example term are not ported (the Julia counterpart of the C
  template is `examples/MyTermPackage/`), and a **"Validation against R"**
  section describing the golden fixture (also in `docs/src/guide/validation.md`).
- **CI derives its sibling clone list from `[sources]`** (`CI.yml` from
  `Project.toml`, `Documentation.yml` from `docs/Project.toml` minus the self
  entry) instead of a hand-written `for pkg in Networks ERGM`; a testset runs
  the same pipeline and asserts it equals `keys([sources])`. The ubuntu/`1`
  cell now also instantiates the benchmark environment and runs
  `regression_tests.jl`, and runs the template package's own `Pkg.test()`
  (`examples/MyTermPackage`), whose `[compat]` now pins `ERGM`,
  `ERGMUserterms` and `Networks` at `0.2` so copied packages pin theirs.

### Changed

- Documentation uses the default Documenter themes, with a new package-specific
  SVG icon and browser favicon in the official Julia logo colors.
- Bundled example terms read attributes per element
  (`get_edge_attribute(net, attr, i, j)` / `get_vertex_attribute(net, attr,
  v)`, which canonicalise the undirected `(min, max)` key) with a default for
  absent values; `DyadCovTerm` reads its matrix at the `(min, max)` entry on
  undirected networks.
- `consistency_check` visits each **unordered** dyad once on an undirected
  network: `exhaustive=true` enumerates `i < j` (the rule `_term_fingerprint`
  already used) and the random budget is `min(100, number of dyads)` with
  draws canonicalised to `(min, max)`. It used to enumerate every ordered
  pair regardless of directedness — 40 `change_stat` calls for the 10 dyads
  of a 5-vertex undirected network, twice the brute-force work — and sized
  the random budget by `n(n-1)`. The verdict is unchanged.
- `validate_traits`' final step reports `ERGMModel construction rejected the
  term or the network: …` (was `… rejected the term on this network: …`);
  the network-only refusals are now raised up front, so what reaches this
  step is a declaration problem, but the wording no longer blames the term
  for what may be the network. A ✓ line is printed only for a check that
  actually ran: the dyad-independence check says so when the network has
  no other dyad to toggle (a 2-vertex undirected network).
- Documentation (term interface, templates, validation) rewritten around
  the add-direction convention, including an explicit migration warning.
- **The validators are deterministic given `rng`.** Every draw comes from
  the `rng` keyword (default `Random.default_rng()`, the task-local RNG, so
  `Random.seed!` still governs unseeded calls); passing `rng=Xoshiro(seed)`
  makes a call independent of global state, and every failure reports the
  literal that replays it. The random-mode `consistency_check` draws up to
  `min(100, n(n-1))` distinct off-diagonal dyads (was `min(100, n²)` pairs
  including self-pairs, which it then skipped).
- `docs/src/index.md` no longer claims "equivalent functionality" to R's
  `ergm.userterms`: the package is its Julia *counterpart in role*, not a
  port of its C contents (which the "Not implemented" section right below
  already said). `getting_started.md`'s "Step 4" and "Complete Example"
  blocks, which re-list their `using` lines, now include the `using Random`
  their `Xoshiro(…)` needs — pasted standalone they raised `UndefVarError:
  Xoshiro`; a docs-gate check now requires every block that lists `using`
  lines and calls `Xoshiro(` to `using Random` itself.
- The edge-attribute snapshot/restore around harness toggles goes through
  Networks.jl's public attribute API (`set_edge_attribute!(…;
  require_edge=false)`) instead of reaching into `net.edge_attrs`.
- The docs describe ERGM 0.2's term contract: a *declared* vertex attribute
  must have a value on every vertex or `ERGMModel` throws (ERGM's message is
  quoted; the check is the public `ERGM._validate_formula`), an undeclared
  read is a `validate_traits` failure, and a `get(attrs, v, default)`
  fallback is reached only by raw `compute` calls (term_interface.md
  "Missing Attribute Values (NA)", templates.md "Using Attributes");
  `Networks.supports_missing === ERGM.supports_missing`, shown with the same
  `ObservedEdges` idiom as ERGM's own docstring (term_interface.md, README,
  template); `Degree`/`Kstar`/`GWDegree` named as ERGM's built-in
  `requires_undirected` examples with the message that names
  `OStar`/`IStar`/`GWODegree`/`GWIDegree`; `name(term, net)` mentioned as
  ERGM's direction-aware label (a user term needs only `name(term)`); the
  workflow table's `ergm()` reads `fit_ergm`/`ergm` (`ergm === fit_ergm`).
  The templates page's local copies of the bundled terms now declare the
  same traits as the shipped ones (`is_dyad_dependent = false`;
  `InteractionTerm`'s two required attributes) and every adapted template
  is followed by an asserted `validate_term`.
- README, `docs/src/index.md` and `getting_started.md` teach one import
  idiom, `import ERGM: name, compute, change_stat` (README and
  `term_interface.md` said `import ERGMUserterms:`; the same functions, and
  the prose now says so once).
- `benchmark_term` is documented as a timer (it times `n_iter` calls; it does
  not profile), with a parameter table; `validation.md`'s "Early Termination"
  matches the code (the interface checks all run; the consistency and trait
  stages are gated on everything before them; `verbose=false` changes only
  the logging).
- The test suite uses ERGM.jl's public trait names (`requires_directed`,
  `required_vertex_attributes`, `has_dyad_dependent`); the sole private
  alias still pinned is `ERGM._requires_directed === ERGM.requires_directed`,
  because TERGM.jl declares methods on it.

### Fixed

- Documentation output comments that the shown code did not produce:
  `validation.md`'s worked example said `compute() returns: 30.0` (the
  value is 26.0) and its `test_term` sample output `Works on complete
  network: 4200.0` (990.0 for `ExampleTerm` on the 10-vertex complete
  digraph; the unseeded random-network value is now shown as `…`);
  `getting_started.md` said `4.0` for `SharedNeighborTerm` on the Step 2
  network (3.0) and `12.0` for the seeded Step 3 run (6.0).
- Three guide pages still taught the whole-Dict
  `get(get_vertex_attribute(net, attr), v, 0.0)` read inside `change_stat`
  that the same pages say never to use: `benchmarking.md` labelled it
  `FAST: O(1)` (it allocates an empty `Dict` per call and infers `Any`; the
  page now shows it as the slow variant and the per-vertex getter as fast),
  and the executed `MyAttrSum` (`term_interface.md`) and `DiffInteraction`
  (`templates.md`) examples now read per element with a `::Float64`
  assertion, as the shipped `InteractionTerm` does. The README's
  `WeightedEdges` sketch (`sum(get_edge_attribute(net, t.attr))` — a
  `MethodError` if pasted, and it ignored `default`) and "Using Attributes"
  pattern (`attrs[src(e)]`, a `KeyError` on any vertex without a value)
  are rewritten with the per-element getters. The template package README's
  "Checking your work" block referenced a `net` no block on the page built
  (`UndefVarError` when pasted); it now builds the test network first.
- `WeightedEdges.compute` sums over the network's actual edges with a
  default weight, so `compute` and `change_stat` stay consistent when the
  sampler adds edges without a weight attribute.
- The validation harness snapshots and restores dyad edge attributes around
  its toggle test, preventing spurious failures for attribute-reading
  terms.
- `test_term`'s `n_tests` keyword was documented but ignored (the network was
  validated with the hard-wired 10 + 5 dyads); it is now passed to
  `validate_term`.
- `consistency_check`'s random mode drew `(rand(1:n), rand(1:n))` pairs,
  so `i == j` self-pairs consumed part of its budget (silently skipped) and
  the draws were not uniform over dyads; it now draws off-diagonal dyads
  directly.
- `validate_term`'s `change_stat` type check reported a failing dyad with no
  way to reproduce it; every failure now carries the replay literal.
- **Docs `using …, Network`** (the pre-rename module name) on four guide
  pages — `benchmarking.md`, `templates.md`, `term_interface.md`,
  `validation.md` (panel 2026-09, item 2) — now reads `Networks`, so the
  worked examples run as written under `tools/check_snippets.jl`; a grep
  guard in the test suite keeps it out.
- **Quick Start examples that could not run as written.**
  `docs/src/getting_started.md`'s walkthrough term iterated `edges(net)` with
  `src`/`dst`, which Networks.jl does not re-export (`UndefVarError` inside
  `compute`, reported by `validate_term` as `compute() failed` — it now
  imports `Graphs: src, dst`); `docs/src/index.md`'s Quick Start defined
  `name`/`compute`/`change_stat` after a bare `using`, so they were *local*
  functions ERGM.jl never called and `validate_term` warned
  `compute() not implemented` while the page "passed" (it now does
  `import ERGM: name, compute, change_stat`, and every such claim is followed
  by `@assert valid`).
- Stale prose: `term_interface.md` showed `name(t::GWDegree) =
  "gwdegree.$(t.decay)"` (ERGM's label is `"gwdeg.fixed.0.5"`, R's); its
  "Missing Attributes" section and `templates.md`'s "Using Attributes" told
  authors to decide for themselves what a missing vertex value means, when
  ERGM 0.2 refuses a declared attribute that is not set on every vertex
  (statnet refuses NA).

## [0.1.0] - 2026-02-09

Initial release: `@ergm_term` scaffolding, example terms, and the term
validation/benchmark harness.
