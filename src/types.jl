
using EnumX

@enumx Efield begin
    "Right Circular Polarization" 
    R
    "Left Circular Polarization"
    L
    "Horizontal Linear Polarization"
    X
    "Vertical Linear Polarization"
    Y
end

export Efield




"""
    $(TYPEDEF)

Denotes a general polarization basis, with basis vectors (B1,B2) which are typically
`<: Union{ElectricFieldBasis, Missing}`
"""
struct PolBasis{P<:Efield.T} 
    p1::P
    p2::P
end


"""
    CirBasis <: PolBasis

Measurement uses the circular polarization basis, which is typically used for circular
feed interferometers. The order is RPol, LPol.
"""
const CirBasis = PolBasis(Efield.R, Efield.L)

const LinBasis = PolBasis(Efield.X, Efield.Y)


"""
    $(TYPEDEF)

Static vector that holds the stokes parameters of a polarized
complex visibility

To convert between a `StokesParams` and `CoherencyMatrix` use the `convert`
function

```julia
convert(::CoherencyMatrix, StokesVector(1.0, 0.1, 0.1, 0.4))
```
"""
struct StokesParams{T} <: FieldVector{4,T}
    I::T
    Q::T
    U::T
    V::T
end

StokesParams(::SArray) = throw(ArgumentError("argument does not have a basis please wrap it in a `CoherencyMatrix`"))


StaticArraysCore.similar_type(::Type{<:StokesParams}, ::Type{T}, s::Size{(4,)}) where {T} = StokesParams{T}


"""
    $(TYPEDEF)

Coherency matrix for a single baseline with bases `B1` and `B2`. The two bases correspond
to the type of feeds used for each telescope and should be subtypes of `PolBasis`. To see which
bases are implemented type `subtypes(Rimes.PolBasis)` in the REPL.

For a circular basis the layout of the coherency matrix is
```
RR* RL*
LR* RR*
```
which can be constructed using
```julia-repl
c = CoherencyMatrix(RR, LR, RL, LL, CirBasis())
```

For a linear basis the layout of the coherency matrix is
```
XX* XY*
YX* YY*
```
which can be constructed using
```julia-repl
c = CoherencyMatrix(XX, YX, XY, YY, CirBasis())
```

For a mixed (e.g., circular and linear basis) the layout of the coherency matrix is
```
RX* RY*
LX* LY*
```

or e.g., linear and circular the layout of the coherency matrix is
```
XR* XL*
YR* YL*
```

These coherency matrices can be constructed using:
```julia-repl
# Circular and linear feeds i.e., |R><X|
c = CoherencyMatrix(RX, LX, RY, LY, LinBasis(), CirBasis())
# Linear and circular feeds i.e., |X><R|
c = CoherencyMatrix(XR, YR, XL, YL, LinBasis(), CirBasis())
```

"""
struct CoherencyMatrix{T, B<:NTuple{2, <:PolBasis}} <: StaticArraysCore.FieldMatrix{2,2,T}
    e11::T
    e21::T
    e12::T
    e22::T
    basis::B
end

StaticArraysCore.similar_type(::Type{CoherencyMatrix{T, B}}, ::Type{T2}, s::Size{(2,2)}) where {T, T2, B} = CoherencyMatrix{T2, B}


"""
    CoherencyMatrix(e11, e21, e12, e22, basis::NTuple{2, PolBasis})

Constructs the coherency matrix with components
   e11 e12
   e21 e22
relative to the tensor product basis, `|basis[1]><basis[2]|`. Note that basis[1] and basis[2]
could be different.

For instance
```julia
c = Coherency(1.0, 0.0, 0.0, 1.0, CirBasis(), LinBasis())
```
elements correspond to
    RX* RY*
    LX* LY*
"""
@inline function CoherencyMatrix(e11::Number, e21::Number, e12::Number, e22::Number, basis::NTuple{2,<:PolBasis})
    T = promote_type(typeof(e11), typeof(e12), typeof(e21), typeof(e22))
    return CoherencyMatrix{T, typeof(basis)}(T(e11), T(e21), T(e12), T(e22), basis)
end

"""
    CoherencyMatrix(e11, e21, e12, e22, basis::PolBasis)

Constructs the coherency matrix with components
   e11 e12
   e21 e22
relative to the tensor product basis, `basis` given by `|basis><basis|`.

For instance
```julia
c = Coherency(1.0, 0.0, 0.0, 1.0, CirBasis())
```
elements correspond to
    RR* RL*
    LR* LL*

"""
@inline function CoherencyMatrix(e11::Number, e21::Number, e12::Number, e22::Number, basis::PolBasis)
    return CoherencyMatrix(e11, e21, e12, e22, (basis, basis))
end

"""
    CoherencyMatrix(e11, e21, e12, e22, basis1::PolBasis basis2::PolBasis)

Constructs the coherency matrix with components
   e11 e12
   e21 e22
relative to the tensor product basis, `basis` given by `|basis1><basis2|`.

For instance
```julia
c = Coherency(1.0, 0.0, 0.0, 1.0, CirBasis(), LinBasis())
```
elements correspond to
    RX* RY*
    LX* LY*

"""
@inline function CoherencyMatrix(e11::Number, e21::Number, e12::Number, e22::Number, basis1::PolBasis, basis2::PolBasis)
    return CoherencyMatrix(e11, e21, e12, e22, (basis1, basis2))
end

@inline function CoherencyMatrix(mat::AbstractMatrix, basis1::PolBasis, basis2::PolBasis)
    return CoherencyMatrix(mat[1], mat[2], mat[3], mat[4], basis1, basis2)
end

@inline function CoherencyMatrix(mat::AbstractMatrix, basis::PolBasis)
    return CoherencyMatrix(mat[1], mat[2], mat[3], mat[4], basis)
end



"""
    CoherencyMatrix(s::StokesParams, basis1::PolBasis)
    CoherencyMatrix(s::StokesParams, basis1::PolBasis, basis2::PolBasis)
    CoherencyMatrix(s::StokesParams, basis1::PolBasis, basis2::PolBasis, refbasis=CirBasis())

Constructs the coherency matrix from the set of stokes parameters `s`.
This is specialized on `basis1` and `basis2` which form the tensor product basis
`|basis1><basis2|`, or if a single basis is given then by `|basis><basis|`.

For example
```julia
CoherencyMatrix(s, CircBasis())
```
will give the coherency matrix
```
   I+V   Q+iU
   Q-iU  I-V
```

while
```julia
CoherencyMatrix(s, LinBasis())
```
will give

```
    I+Q   U+iV
    U-iV  I-Q
```

# Notes

Internally this function first converts to a reference basis and then the final basis.
You can select the reference basis used with the optional argument refbasis. By default
we use the circular basis as our reference. Note that this is only important for mixed bases,
e.g., if `basis1` and `basis2` are different. If `basis1==basis2` then the reference basis
is never used.
"""
@inline function CoherencyMatrix(s::StokesParams, b1::PolBasis, b2::PolBasis, refbasis=CirBasis)
    t1 = basis_transform(refbasis=>b1)
    # Flip because these are the dual elements
    t2 = basis_transform(b2=>refbasis)
    # Use circular basis as a reference
    c_cir = CoherencyMatrix(s, refbasis)
    return CoherencyMatrix(t1*c_cir*t2, b1, b2)
end

@inline function CoherencyMatrix(s::StokesParams, b::PolBasis)
    if b == CirBasis || b == PolBasis(Efield.L, Efield.R)
        (;I,Q,U,V) = s
        RR = complex((I + V))
        LR = (Q - 1im*U)
        RL = (Q + 1im*U)
        LL = complex((I - V))
        b == CirBasis && return CoherencyMatrix(RR, LR, RL, LL, (b,b))
        return CoherencyMatrix(LL, RL, LR, RR, (b,b))
    elseif b == LinBasis || b == PolBasis(YPol, XPol)
        (;I,Q,U,V) = s
        XX = complex((I + Q))
        YX = (U - 1im*V)
        XY = (U + 1im*V)
        YY = complex((I - Q))
        b == CirBasis && return CoherencyMatrix(XX, YX, XY, YY, (b,b))
        return CoherencyMatrix(YY, XY, YX, XX, (b,b))
    else
        throw(ArgumentError("Unsupported basis $b"))
    end
end


@inline function StokesParams(c::CoherencyMatrix)
    if c.basis[1] == c.basis[2] == CirBasis
        I = (c.e11 + c.e22)/2
        Q = (c.e21 + c.e12)/2
        U = 1im*(c.e21 - c.e12)/2
        V = (c.e11 - c.e22)/2
        return StokesParams(I, Q, U, V)
    else
        return _stokesparams(c)
    end
end

@inline function _stokesparams(c::CoherencyMatrix)
    b = c.basis
    t1 = basis_transform(b[1]=>CirBasis)
    # Flip because these are the dual elements
    t2 = basis_transform(CirBasis=>b[2])
    c_cir = CoherencyMatrix(t1*c*t2, CirBasis)
    return StokesParams(c_cir)
end
