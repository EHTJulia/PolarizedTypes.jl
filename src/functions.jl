"""
    $(SIGNATURES)

Computes `linearpol` from a set of stokes parameters `s`.
"""
function linearpol(s::StokesParams{T}) where {T}
    Tr = real(T)
    im = complex(zero(Tr), one(Tr))
    return s.Q + im * s.U
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
    a = (p + abs(l)) / 2
    b = (p - abs(l)) / 2
    ev = evpa(s)
    sn = sign(s.V)
    return (; a, b, evpa = ev, sn)
end


"""
    $(SIGNATURES)

Returns the (Q/I, U/I, V/I) fractional polarization vector as a 3-element static vector.
"""
fracpolarization(s::StokesParams) = polarization(s) / s.I

"""
    mpol(m::StokesParameters{<:Real})


Compute the complex fractional linear polarization of a Stokes Parameter `m`
"""
mpol(m::StokesParams) = linearpol(m) / m.I

"""
    evpa(m::Union{StokesParams, CoherencyMatrix})
Compute the evpa of a stokes vect or cohereny matrix.
"""
evpa(m::StokesParams) = angle(linearpol(m)) / 2


# Now define the functions for Coherency matrix
for f in (:linearpol, :polarization, :fracpolarization, :mpol, :evpa)
    @eval begin
        $(f)(c::CoherencyMatrix) = $(f)(StokesParams(c))
    end
end
