"""
    $(SIGNATURES)

Computes `linearpol` from a set of stokes parameters `s`.
"""
function linearpol(s::StokesParams)
    return s.Q + 1im*s.U
end

"""
    $(SIGNATURES)

Returns the (Q, U, V) polarization vector as a 3-element static vector.
"""
function polarization(s::StokesParams)
    return SVector(s.Q, s.U, s.V)
end

"""
    $(SIGNATURES)

Returns the polarization ellipse of the Stokes parameters `s`. The results
is a named tuple with elements
  - `a`   : The semi-major axis of the polarization ellipse
  - `b`   : The semi-minor axis of the polarization ellipse
  - `evpa`: The electric vector position angle of `s` or the PA of the ellipse.
  - `sn`  : The sign of the Stokes `V`.


## Notes

The semi-major and semi-minor axes are defined as
```
    a = 1/2(Iₚ + |L|)
    b = 1/2(Iₚ - |L|)
```
where `Iₚ = √(Q² + U² + V²)` and `|L| = √(Q² + U²)`.

In general the area of the ellipse is given by
```
    πab = π/4|V|²
```

For sources with zero linear polarization `a = b` so we
have a circle with radius `|V|`. For purely linear polarization `b = 0` giving a line
with length `|L|`.

"""
function polellipse(s::StokesParams{<:Real})
    l = linearpol(s)
    p = norm(polarization(s))
    a = (p + abs(l))/2
    b = (p - abs(l))/2
    ev = evpa(s)
    sn = sign(s.V)
    return (;a, b, evpa=ev, sn)
end



"""
    $(SIGNATURES)

Returns the (Q/I, U/I, V/I) fractional polarization vector as a 3-element static vector.
"""
fracpolarization(s::StokesParams) = polarization(s)*inv(s.I)

"""
    mpol(m::StokesParameters{<:Real})


Compute the complex fractional linear polarization of a Stokes Parameter `m`
"""
mpol(m::StokesParams) = (m.Q + 1im*m.U)/m.I

"""
    m̆(m::Union{StokesParameters{<:Complex}, CoherencyMatrix)

Computes the complex fractional linear polarization of the complex or visibility quantities.
Note that this function can also be called used [`mbreve`](@ref) or [`mpol`](@ref).
"""
m̆(m::StokesParams) = mpol(m)


"""
    $(SIGNATURES)

Computes the complex fractional linear polarization of the complex or visibility quantities.
Note that this function can also be called used [`m̆`](@ref) or [`mpol`](@ref).
"""
mbreve(m::StokesParams) = m̆(m)

"""
    evpa(m::Union{StokesParams, CoherencyMatrix})
Compute the evpa of a stokes vect or cohereny matrix.
"""
evpa(m::StokesParams) = 1/2*atan(m.U, m.Q)
evpa(m::StokesParams{<:Complex}) = 1/2*angle(m.U/m.Q)


# Now define the functions for Coherency matrix
for f in (:linearpol, :polarization, :fracpolarization, :m̆, :mbreve, :evpa)
    @eval begin
        $(f)(c::CoherencyMatrix) = $(f)(StokesParams(c))
    end
end
