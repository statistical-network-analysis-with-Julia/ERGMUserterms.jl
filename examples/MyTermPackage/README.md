# MyTermPackage — a template for third-party ERGM terms

A complete, runnable skeleton of a package that ships a custom ERGM term. Copy
the directory, rename the module, replace the term, keep the structure.

```
MyTermPackage/
├── Project.toml        # deps: ERGM, ERGMUserterms, Networks, Graphs
├── src/
│   └── MyTermPackage.jl   # the term: 3 interface methods + 4 trait declarations
└── test/
    └── runtests.jl        # validate_term + ERGMModel construction + fit (seeded)
```

## What a term owes ERGM.jl

**The interface** (`import ERGM: compute, change_stat, name` — extend the
shared generics, never define local ones):

| method                       | contract                                                          |
|:-----------------------------|:------------------------------------------------------------------|
| `name(t)`                    | lowercase, parameterized, R-ergm style (`"recip_homophily.group"`) |
| `compute(t, net)`            | the full statistic `g(y)`                                          |
| `change_stat(t, net, i, j)`  | **add-direction** `g(y⁺ᵢⱼ) − g(y⁻ᵢⱼ)`, *independent of the dyad's current state* |

The toggle-direction idiom `has_edge(net, i, j) ? -Δ : Δ` is wrong and
`validate_term` rejects it.

**The traits** — ERGM.jl's public term-trait protocol. They are how your term
tells ERGM what data it needs and which models it belongs in; ERGM acts on them
when an `ERGMModel` is built, so a third-party term is validated exactly like a
built-in one:

```julia
ERGM.required_vertex_attributes(t::MyTerm) = (t.attr,)   # default ()
ERGM.required_edge_attributes(t::MyTerm)   = ()          # default ()
ERGM.requires_directed(::MyTerm)           = true        # default false
ERGM.requires_undirected(::MyTerm)         = false       # default false
ERGM.is_dyad_dependent(::MyTerm)           = true        # default true (conservative)
Networks.supports_missing(::MyTerm)        = true        # default false; ≡ ERGM.supports_missing
```

- **Declare an attribute** iff its absence is an error. `get_vertex_attribute`
  returns an empty `Dict` for an unknown attribute, so an undeclared
  attribute-based term silently becomes an all-zero design column and a
  meaningless coefficient; declared, it raises an `ArgumentError` naming the
  attribute. Do *not* declare an attribute the term merely reads *if present*
  (a weight with a default, say) — that would reject networks the term handles
  fine.
- **Declare a direction requirement** iff the statistic is undefined (or a
  duplicate of another term) on the other kind of network.
- **Declare `is_dyad_dependent = false`** iff the change statistic never reads
  another dyad's state. It buys you an exact MPLE and no pseudo-likelihood
  caveat — claim it only if it is true.
- **Declare `supports_missing = true`** iff the statistic *consults the mask*
  (`is_missing_dyad`) and is therefore invariant to the face value of an
  unobserved dyad. The default `false` is the honest answer for terms that
  simply count edges as stored. `Networks.supports_missing` and
  `ERGM.supports_missing` are the same function, so either spelling extends
  it.
- **A declared vertex attribute must be complete.** `ERGMModel` refuses a
  network on which it is set on only some vertices (`… needs vertex attribute
  :group on every vertex … statnet refuses NA`), so the `get(vals, v,
  nothing)` fallback in `_counts` is reached only by raw `compute` calls.

## Checking your work

On the template's own test network (the one `test/runtests.jl` builds —
8 vertices, a `:group` attribute on every vertex, seven arcs):

```julia
using ERGM, ERGMUserterms, Networks, MyTermPackage, Random
net = network(8; directed=true)
set_vertex_attribute!(net, :group, Dict(v => isodd(v) ? "a" : "b" for v in 1:8))
for (i, j) in ((1, 3), (3, 1), (2, 4), (4, 2), (1, 2), (5, 7), (6, 8))
    add_edge!(net, i, j)
end
validate_term(ReciprocatedHomophily(:group), net; rng=Xoshiro(1))   # interface AND traits
```

The harness draws every random dyad from `rng` (default
`Random.default_rng()`); pass a seeded `Xoshiro` in your tests so a red run
is the same run every time. A failure message ends with
`reproduce with rng=Xoshiro(0x…, …) passed to validate_term` — hand it back
to replay that exact run.

`validate_term` brute-force-toggles dyads to verify the add-direction
convention, perturbs each attribute to check that the declared ones are the
ones the term actually reads, builds the direction twin of the network to check
that a direction requirement is enforced, toggles other dyads to check a
dyad-independence claim, and flips the face value of masked dyads to check a
`supports_missing` claim. Then run your term through a real fit — with an
`rng`, always (every fit in the ecosystem takes one; MPLE happens to be
deterministic, `method=:mcmle` is not):

```julia
using ERGM, Random
fit_ergm(net, [Edges(), ReciprocatedHomophily(:group)]; rng=Xoshiro(1))
```

## Running the template's tests

```bash
cd examples/MyTermPackage
julia --project -e 'using Pkg; Pkg.test()'
```

(The same term and assertions are exercised from ERGMUserterms.jl's own test
suite, so the template cannot rot.)
