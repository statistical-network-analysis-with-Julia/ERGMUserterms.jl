# ERGMUserterms.jl

*Custom ERGM Term Development for Julia*

A Julia package providing templates, utilities, and validation tools for developing custom ERGM terms.

## Overview

ERGMUserterms.jl is a toolkit for creating, validating, and benchmarking custom terms for Exponential Random Graph Models (ERGMs). It provides a structured development workflow with automatic consistency checking, ensuring that user-defined terms correctly implement the ERGM term interface before they are used in model estimation.

ERGMUserterms.jl is a port of [ergm.userterms](https://github.com/statnet/ergm.userterms) from the StatNet collection, providing equivalent functionality in Julia.

### What is a Custom ERGM Term?

An ERGM term is a network statistic that captures a structural feature:

```text
compute(term, net) → network statistic value
change_stat(term, net, i, j) → add-direction change: g(y⁺ij) − g(y⁻ij)
```

Custom terms let you model network features beyond the built-in terms in ERGM.jl, such as:

- Domain-specific structural patterns
- Weighted or attributed network effects
- Novel configurations involving node or edge covariates
- Interaction effects between multiple attributes

### Key Concepts

| Concept | Description |
|---------|-------------|
| **ERGM Term** | A network statistic measuring a structural feature |
| **compute()** | Calculates the full statistic value for a network |
| **change_stat()** | The add-direction change statistic g(y⁺ij) − g(y⁻ij), state-independent |
| **Consistency** | `change_stat(net, i, j) == compute(net') - compute(net)` |
| **Validation** | Automated checks that a term is correctly implemented |

### Use Cases

ERGMUserterms.jl is designed for:

- **Methodologists**: Developing novel ERGM terms for new research questions
- **Applied researchers**: Creating domain-specific terms not available in ERGM.jl
- **Package developers**: Building ERGM extensions with validated term implementations
- **Students**: Learning the ERGM term interface through templates and examples

## Features

- **Term macro**: `@ergm_term` for defining custom terms with automatic validation
- **Validation framework**: Comprehensive checks for `name()`, `compute()`, and `change_stat()` correctness
- **Consistency checking**: Verify that `change_stat()` matches `compute()` differences
- **Performance profiling**: Benchmark `compute()` and `change_stat()` to ensure efficiency
- **Example terms**: Five complete, working examples to copy and modify

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/Statistical-network-analysis-with-Julia/ERGMUserterms.jl")
```

Or for development:

```julia
using Pkg
Pkg.develop(path="/path/to/ERGMUserterms.jl")
```

## Quick Start

```julia
using ERGM
using ERGMUserterms
using Network

# Define a custom term
struct MyTerm <: AbstractUserTerm
    param::Float64
end

name(t::MyTerm) = "myterm.$(t.param)"

function compute(t::MyTerm, net)
    return Float64(ne(net)) * t.param
end

function change_stat(t::MyTerm, net, i::Int, j::Int)
    # Add-direction: statistic with edge (i,j) present minus absent,
    # regardless of the dyad's current state
    return t.param
end

# Create a test network and validate
net = network(20; directed=true)
validate_term(MyTerm(2.0), net)
```

## Term Development Workflow

| Step | Function | Description |
|------|----------|-------------|
| 1. Define | `struct MyTerm <: AbstractUserTerm` | Create term struct with fields |
| 2. Implement | `name()`, `compute()`, `change_stat()` | Implement the three required methods |
| 3. Validate | [`validate_term`](@ref) | Check all methods work correctly |
| 4. Test | [`test_term`](@ref) | Run comprehensive tests on random networks |
| 5. Benchmark | [`benchmark_term`](@ref) | Profile performance of compute vs change_stat |
| 6. Use | Pass to ERGM.jl | Use in `ergm()` model specification |

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
