# ERGMUserterms.jl

Develop a network statistic for ERGM.jl and verify that its change statistic agrees with full recomputation. ERGMUserterms.jl provides a Julia term interface, validation tools, benchmarks, and a copyable package template for extending binary ERGMs.

| First analysis | Learn the model or data | Reference and detail |
|:--|:--|:--|
| [Implement a worked term](getting_started.md) | [Validate the implementation](guide/validation.md) | [Browse the term interface](api/types.md) |

!!! note "Supported scope"

    The validator checks computational consistency on the networks it exercises. It does not establish identifiability, model fit, or scientific validity. A custom term must declare its dependence and attribute requirements accurately; missing-data support is opt-in.

## Installation

```@raw html
<p>Use Julia <strong>1.12 or newer</strong> and the <a href="/getting-started/">shared workspace installation guide</a>. These development packages are not yet registered; the guide prepares the required sibling checkouts and a Julia environment for the examples.</p>
```

## Quick Start

Implement a scaled edge count. Import the shared function names before adding methods:

```julia
using NetworkCore, ERGM, ERGMUserterms, Random
import ERGM: name, compute, change_stat

struct ScaledEdgesPreview <: AbstractUserTerm
    scale::Float64
end
name(t::ScaledEdgesPreview) = "scaled_edges.$(t.scale)"
compute(t::ScaledEdgesPreview, net::Network) = t.scale * ne(net)
change_stat(t::ScaledEdgesPreview, net::Network, i::Int, j::Int) = t.scale
ERGM.is_dyad_dependent(::ScaledEdgesPreview) = false

net = load_dataset(:florentine_marriage)
@assert validate_term(ScaledEdgesPreview(2.0), net; rng=Xoshiro(1))
```

`change_stat` always returns the addition-direction difference, including when the edge already exists; the sampler handles the removal sign. This statistic is proportional to `Edges()`, so using both in one model would create redundant parameters.

## Not implemented (vs R ergm.userterms)

What is ported is the *role* of `ergm.userterms` — a validated, copyable
starting point for third-party terms — not its contents. R's package is a C
template: its `changestats.users.c` skeleton and its worked example term
`mindegree` are **not** ported, and nothing shipped here has them as a
counterpart. A Julia term is plain Julia (three methods plus trait
declarations), so there is no C skeleton to fill in; the Julia counterpart of
the C template is the package template `examples/MyTermPackage/`, and the
bundled example terms are instead pinned against the ergm terms they are
identical to (see [Validation against R](@ref)).

Limits of what the validation harness establishes:

- **A pass is a statement about the network you validated on.** By default
  every dyad of that network is checked (and the path to it from the empty
  network), but a change statistic that is wrong only in a configuration the
  network does not contain is not seen. Validate on networks that contain
  the configurations your term counts, and use [`test_term`](@ref) for random ones.
- **Above 5000 dyads the per-dyad checks sample** (stratified by dyad state,
  5000 dyads), and the dyad-independence check toggles every pair of dyads
  only up to about 500 dyads (22 vertices directed, 32 undirected); beyond
  that it toggles the reverse arc and the dyads sharing an endpoint first and
  a random remainder. The log line states how many dyads were checked.
- **Edge attributes are not preserved by ERGM.jl's samplers.** A term that
  reads one must hold it as a matrix (as `WeightedEdges` and `DyadCovTerm`
  do). `validate_traits` fails a term that reads it live, and ERGM.jl's
  samplers (`mcmle`, `simulate_ergm`, `gof`, the MPLE bootstrap) refuse one
  with an `ArgumentError`; the MPLE, which never toggles a tie, accepts it.
  ERGM.jl does trust a term's `is_dyad_dependent` declaration beyond a
  probe at model construction, which is why the validator tests that claim.
- The harness checks computation, not statistics: it does not establish
  identifiability, non-degeneracy or model fit.

## Term Development Workflow

| Step | Function | Description |
|------|----------|-------------|
| 1. Define | `struct MyTerm <: AbstractUserTerm` | Create term struct with fields |
| 2. Implement | `name()`, `compute()`, `change_stat()` | Implement the three required methods |
| 3. Validate | [`validate_term`](@ref) | Check all methods work correctly |
| 4. Test | [`test_term`](@ref) | Run comprehensive tests on random networks (`directed=false` / `vertex_attributes=` for terms that need them) |
| 5. Benchmark | [`benchmark_term`](@ref) | Time compute vs change_stat |
| 6. Use | Pass to ERGM.jl | Put it in the term list of `fit_ergm`/`ergm` (`ergm === fit_ergm`, a `const` alias) |

## Documentation

```@contents
Pages = [
    "getting_started.md",
    "guide/term_interface.md",
    "guide/templates.md",
    "guide/validation.md",
    "guide/benchmarking.md",
    "api/types.md",
    "api/validation.md",
    "api/utilities.md",
]
Depth = 2
```

## Theoretical Background

### The ERGM Term Interface

In an Exponential Random Graph Model, the probability of a network $\mathbf{Y}$ is:

$$P(\mathbf{Y} = \mathbf{y}) = \frac{1}{\kappa(\boldsymbol{\theta})} \exp\left(\sum_k \theta_k g_k(\mathbf{y})\right)$$

Where:

- $g_k(\mathbf{y})$ are network statistics (terms) computed by `compute()`
- $\theta_k$ are parameters to be estimated
- $\kappa(\boldsymbol{\theta})$ is the normalizing constant

### Change Statistics in MCMC

ERGM estimation relies on MCMC simulation, where edges are toggled one at a time. The change statistic for dyad $(i,j)$ is the add-direction difference

$$\Delta g_k(i,j) = g_k(\mathbf{y}^{+}_{ij}) - g_k(\mathbf{y}^{-}_{ij})$$

with $\mathbf{y}^{+}_{ij}$/$\mathbf{y}^{-}_{ij}$ the network with the edge forced present/absent. `change_stat()` must return exactly this quantity — independent of the dyad's current state (the sampler negates it for removals). ERGMUserterms.jl validates both the value and its state-independence.

## References

1. Hunter, D.R., Handcock, M.S., Butts, C.T., Goodreau, S.M., Morris, M. (2008). ergm: A package to fit, simulate and diagnose exponential-family models for networks. *Journal of Statistical Software*, 24(3), 1-29.

2. Morris, M., Handcock, M.S., Hunter, D.R. (2008). Specification of exponential-family random graph models: Terms and computational aspects. *Journal of Statistical Software*, 24(4), 1-24.

3. Hunter, D.R. (2007). Curved exponential family models for social networks. *Social Networks*, 29(2), 216-230.

4. Robins, G., Pattison, P., Kalish, Y., Lusher, D. (2007). An introduction to exponential random graph (p*) models for social networks. *Social Networks*, 29(2), 173-191.

## Citation

If you use ERGMUserterms.jl in your work, please cite it using the entry in
[`CITATION.bib`](https://github.com/statistical-network-analysis-with-Julia/ERGMUserterms.jl/blob/main/CITATION.bib).
Please also cite the R packages it is a counterpart of and their methods
papers — `ergm.userterms` (Hunter, Goodreau & Handcock 2013) and `ergm`
(Hunter et al. 2008; Krivitsky et al. 2023); the ecosystem's
[How to cite](https://statistical-network-analysis-with-julia.github.io/citing/)
page lists the references with DOIs.

```biblatex
@misc{SNWJERGMUsertermsJL,
  author = {Santoni, Simone},
  title = {ERGMUserterms.jl: Custom ERGM Term Development for Julia},
  year = {2026},
  url = {https://github.com/statistical-network-analysis-with-Julia/ERGMUserterms.jl},
  note = {Homepage: https://statistical-network-analysis-with-Julia.github.io/ERGMUserterms.jl; GitHub: https://github.com/statistical-network-analysis-with-Julia}
}
```
