module PolarizedTypes

using DocStringExtensions
using PrecompileTools
using StaticArrays
using LinearAlgebra

# Write your package code here.
export StokesParams, CoherencyMatrix, RPol, LPol, XPol, YPol,
    PolBasis, CirBasis, LinBasis,
    basis_components, basis_transform, coherencymatrix, stokesparams,
    linearpol, mbreve, m̆, evpa, polarization, fracpolarization, mpol,
    polellipse

include("types.jl")
include("basis_transforms.jl")
include("functions.jl")

@setup_workload begin
    @compile_workload begin
        # Now polarization stuff
        s64 = StokesParams(1.0, 0.5, 0.5, 0.5)
        c164 = CoherencyMatrix(s64, CirBasis(), CirBasis())
        c264 = CoherencyMatrix(s64, LinBasis(), LinBasis())
        c364 = CoherencyMatrix(s64, CirBasis(), LinBasis())
        c464 = CoherencyMatrix(s64, LinBasis(), CirBasis())

        StokesParams(c164)
        StokesParams(c264)
        StokesParams(c364)
        StokesParams(c464)

        s32 = StokesParams(1.0f0 + 0im, 5.0f-1 + 0im, 5.0f-1, 5.0f-1)
        c132 = CoherencyMatrix(s32, CirBasis(), CirBasis())
        c232 = CoherencyMatrix(s32, LinBasis(), LinBasis())
        c332 = CoherencyMatrix(s32, CirBasis(), LinBasis())
        c432 = CoherencyMatrix(s32, LinBasis(), CirBasis())

        StokesParams(c132)
        StokesParams(c232)
        StokesParams(c332)
        StokesParams(c432)
    end
end

end
