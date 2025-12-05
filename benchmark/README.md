# Benchmark & Performance

This folder contains scripts and tools to evaluate the performance of **EXOSPHID**.

---

## General Benchmark

The main function for benchmarking is:

```julia
using EXOSPHID

photobenchmark(dt, solar_activity, parent_name)
```

### Function Inputs

- `dt`  
  Time window for photon interaction. Choose sufficiently high to obtain eprformance for instances where photoreactions are actually hapening. 

- `solar_activity`  
  Solar activity level, from `0` (Quiet Sun) to `1` (Active Sun).  For validation against literature values, always use Quiet Sun (`0`).

- `parent_name`  
  Type of parent molecule. Must be one of the **EXOSPHID species**, e.g., `"H2O"`, `"OH"`, `"H2"`. 

### Example

```julia
using EXOSPHID
include(".../benchmark/benchmark.jl")
# Benchmark H2O performance under quiet Sun conditions with 1e10 s time window
photobenchmark(1e10, 0, "H2O")
```

---

## Test Allocations

Use Julia's [`Profile`]() package to track numerical allocations.

### Example

```julia
using EXOSPHID
include(".../benchmark/benchmark.jl")
allocations(1e10, 0, "H2O")
```

---

For more information on performance and benchmarking results see [WIKI: Benchmark]().