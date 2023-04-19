# LocalFennelPartitioning

<!-- [![Build Status](https://github.com/TrainOfCode/LocalFennelPartitioning.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/TrainOfCode/LocalFennelPartitioning.jl/actions/workflows/CI.yml?query=branch%3Amain) -->

This package contains the function
```
local_fennel(G::Vector{Vector{Int}}, locations::Matrix{Float64}, num_p::Int, alpha::Float64, gamma::Float64, to_cons::Int)
```

and

```
local_fennel_sim(G::Vector{Vector{Int}}, locations::Matrix{Float64}, num_p::Int)
```

In `local_fennel_sim`; `alpha`, `gamma`, and `to_cons` are calculated internally

The two functions partition G using a modified version very similar to the Fennel Algorithm, producing the following

```
=> (partitions: Vector{Vector{Int}}, lookups:Vector{Int})
  where for node i ∈ Vertices(G), node i has been assigned to the partition number lookups[i], and likewise i ∈ partitions[lookups[i]]
```
