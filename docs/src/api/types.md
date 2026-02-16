# Types API Reference

This page documents the core types in ERGMUserterms.jl.

## Term Base Types

### AbstractUserTerm

```@docs
AbstractUserTerm
```

## Example Terms

### ExampleTerm

```@docs
ExampleTerm
```

### TemplateTerm

```@docs
TemplateTerm
```

### WeightedEdges

```@docs
WeightedEdges
```

### DyadCovTerm

```@docs
DyadCovTerm
```

### InteractionTerm

```@docs
InteractionTerm
```

## Term Interface

All terms (both user-defined and built-in) must implement these methods.

### name

```@docs
name
```

### compute

```@docs
compute
```

### change_stat

```@docs
change_stat
```
