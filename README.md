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

In `local_fennel_sim`; `alpha`, `gamma`, and `to_cons` are calculated internally, and then it calls local_fennel

Local Fennel is based off of the streaming partitioning algorithm that goes by the name Fennel, I have modified it for my own purposes, and is available here

Both functions produce the following

```
=> (partitions: Vector{Vector{Int}}, lookups:Vector{Int})
  where for node i ∈ Vertices(G), node i has been assigned to the partition number lookups[i], and likewise i ∈ partitions[lookups[i]]
```

Fennel Algorithm: Charalampos Tsourakakis, Christos Gkantsidis, Bozidar Radunovic, and Milan Vojnovic. 2014. FENNEL: streaming graph partitioning for massive scale graphs. In Proceedings of the 7th ACM international conference on Web search and data mining (WSDM '14). Association for Computing Machinery, New York, NY, USA, 333–342. https://doi.org/10.1145/2556195.2556213
