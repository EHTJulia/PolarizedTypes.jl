module ChainRulesCoreExt

using PolarizedTypes
using ChainRulesCore
using ChainRulesCore: ProjectTo, NoTangent

ChainRulesCore.ProjectTo(x::CoherencyMatrix{<: Number}) = ProjectTo{CoherencyMatrix}(; element = ProjectTo(eltype(x)), basis = x.basis)

function (project::ProjectTo{CoherencyMatrix})(dx::AbstractMatrix)
    @assert size(dx) == (2,2) "Issue in Coherency pullback the matrix is not 2x2"
    return CoherencyMatrix(dx, project.basis)
end

function (project::ProjectTo{CoherencyMatrix})(dx::CoherencyMatrix) 
    @assert dx.basis == project.basis "First basis does not match in $(typeof(dx)) and $(project.basis)"
    @assert size(dx) == (2,2) "Issue in Coherency pullback the matrix is not 2x2"
    return dx
end

# function ChainRulesCore.rrule(::Type{<:CoherencyMatrix}, e11, e21, e12, e22, basis::NTuple{2, <:PolBasis})
#     c = CoherencyMatrix(e11, e21, e12, e22, basis)
#     pr = ProjectTo(e11)
#     function _CoherencyMatrix_sep_pullback(Δ)
#         return NoTangent(), pr(Δ[1,1]), pr(Δ[2,1]), pr(Δ[1,2]), pr(Δ[2,2]), NoTangent()
#     end
#     return c, _CoherencyMatrix_sep_pullback
# end

end
