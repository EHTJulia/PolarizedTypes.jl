# PolarizedTypes

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ehtjulia.github.io/PolarizedTypes.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ehtjulia.github.io/PolarizedTypes.jl/dev/)
[![Build Status](https://github.com/ehtjulia/PolarizedTypes.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ehtjulia/PolarizedTypes.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/ehtjulia/PolarizedTypes.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/ehtjulia/PolarizedTypes.jl)
[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)


This defines the basic for polarized types for use in VLBI, including:
  - `StokesParams` for the stokes parameters
  - `CoherencyMatrix` for coherency matrices in arbitrary bases, including a mixed basis.

```julia
using PolarizedTypes

s = StokesParams(1.0, 0.1, 0.1, -0.05)
c = CoherencyMatrix(s, CirBasis(), CirBasis())
l = CoherencyMatrix(s, CirBasis(), LinBasis())
m = CoherencyMatrix(x, LinBasis(), CirBasis())
```
