ChainRulesCore.ProjectTo(x::CoherencyMatrix{B1, B2, <: Number}) where {B1, B2} = ProjectTo{CoherencyMatrix}(; element = ProjectTo(eltype(x)), basis1=B1(), basis2=B2())
function (project::ProjectTo{CoherencyMatrix})(dx::AbstractMatrix)
    @assert size(dx) == (2,2) "Issue in Coherency pullback the matrix is not 2x2"
    return CoherencyMatrix(dx, project.basis1, project.basis2)
end


# function ChainRulesCore.rrule(::Type{<:CoherencyMatrix}, e11, e21, e12, e22, basis::NTuple{2, <:PolBasis})
#     c = CoherencyMatrix(e11, e21, e12, e22, basis)
#     pr = ProjectTo(e11)
#     function _CoherencyMatrix_sep_pullback(Δ)
#         return NoTangent(), pr(Δ[1,1]), pr(Δ[2,1]), pr(Δ[1,2]), pr(Δ[2,2]), NoTangent()
#     end
#     return c, _CoherencyMatrix_sep_pullback
# end

# Needed to ensure everything is constructed nicely
