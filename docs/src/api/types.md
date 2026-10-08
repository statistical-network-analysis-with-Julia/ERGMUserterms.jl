# Types API Reference

This page documents the core types in ERGMUserterms.jl.

## Package

```@docs
ERGMUserterms
```

## Term Base Types

```@docs
AbstractUserTerm
```

## Example Terms

```@docs
ExampleTerm
```

```@docs
TemplateTerm
```

```@docs
WeightedEdges
```

```@docs
DyadCovTerm
```

```@docs
InteractionTerm
```

## Term Interface

All terms (both user-defined and built-in) must implement these methods.
The generics are owned by ERGM.jl and re-exported by ERGMUserterms.jl.

```@docs
name
```

```@docs
compute
```

```@docs
change_stat
```
