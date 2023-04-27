"""
    $(SIGNATURES)

Computes `linearpol` from a set of stokes parameters `s`.
"""
function linearpol(s::StokesParams)
    return s.Q + 1im*s.U
end

"""
    $(SIGNATURES)
Compute the fractional linear polarization of a stokes vector
or coherency matrix
"""
m̆(m::StokesParams{T}) where {T} = (m.Q + 1im*m.U)/(m.I + eps(T))
m̆(m::StokesParams{Complex{T}}) where {T} = (m.Q + 1im*m.U)/(m.I + eps(T))
m̆(m::CoherencyMatrix{CirBasis,CirBasis}) = 2*m.e12/(m.e11+m.e22)
mbreve(m::Union{StokesParams, CoherencyMatrix}) = m̆(m)

"""
    $(SIGNATURES)
Compute the evpa of a stokes vector or cohereny matrix.
"""
evpa(m::StokesParams) = 1/2*atan(m.U, m.Q)
evpa(m::StokesParams{<:Complex}) = 1/2*angle(m.U/m.Q)
evpa(m::CoherencyMatrix{CirBasis, CirBasis}) = angle(m.e12)
