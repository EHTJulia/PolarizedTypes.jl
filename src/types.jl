"""
    $(TYPEDEF)

An abstract type whose subtypes denote a specific electric field basis.
"""
abstract type ElectricFieldBasis end

"""
    $(TYPEDEF)

The right circular electric field basis, i.e. a right-handed circular feed.
"""
struct RPol <: ElectricFieldBasis end

"""
    $(TYPEDEF)

The left circular electric field basis, i.e. a left-handed circular feed.
"""
struct LPol <: ElectricFieldBasis end

"""
    $(TYPEDEF)

The horizontal or X electric feed basis, i.e. the horizontal linear feed.
"""
struct XPol <: ElectricFieldBasis end

"""
    $(TYPEDEF)

The vertical or Y electric feed basis, i.e. the vertical linear feed.
"""
struct YPol <: ElectricFieldBasis end






"""
    $(TYPEDEF)

Denotes a general polarization basis, with basis vectors (B1,B2) which are typically
`<: Union{ElectricFieldBasis, Missing}`
"""
struct PolBasis{B1<:Union{ElectricFieldBasis, Missing}, B2<:Union{ElectricFieldBasis, Missing}} end


"""
    CirBasis <: PolBasis

Measurement uses the circular polarization basis, which is typically used for circular
feed interferometers. The order is RPol, LPol.
"""
const CirBasis = PolBasis{RPol,LPol}

"""
    LinBasis <: PolBasis

Measurement uses the linear polarization basis, which is typically used for linear
feed interferometers. Tge order is XPol, YPol.
"""
const LinBasis = PolBasis{XPol,YPol}

const UnionPolBasis = Union{CirBasis, PolBasis{LPol, RPol}, 
                            LinBasis, PolBasis{YPol, XPol}}

# Horrible hack to automatically promote vectors to use the UnionPolBasis type
# if applicable.
Base.promote_rule(::Type{P}, ::Type{P}) where {P<:PolBasis} = P
Base.promote_rule(::Type{P1}, ::Type{P2}) where {P1 <: UnionPolBasis, P2 <: UnionPolBasis} = UnionPolBasis


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
struct CoherencyMatrix{B1,B2,T} <: StaticArraysCore.FieldMatrix{2,2,T}
    e11::T
    e21::T
    e12::T
    e22::T
end

StaticArraysCore.similar_type(::Type{CoherencyMatrix{B1,B2}}, ::Type{T}, s::Size{(2,2)}) where {B1,B2,T} = CoherencyMatrix{B1,B2,T}


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
@inline function CoherencyMatrix(e11::Number, e21::Number, e12::Number, e22::Number, basis::NTuple{2,PolBasis})
    T = promote_type(typeof(e11), typeof(e12), typeof(e21), typeof(e22))
    return CoherencyMatrix{typeof(basis[1]), typeof(basis[2]),T}(T(e11), T(e21), T(e12), T(e22))
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
@inline function CoherencyMatrix(s::StokesParams, b1::PolBasis, b2::PolBasis, refbasis::Union{LinBasis, CirBasis}=CirBasis())
    t1 = basis_transform(refbasis=>b1)
    # Flip because these are the dual elements
    t2 = basis_transform(b2=>refbasis)
    # Use circular basis as a reference
    c_cir = CoherencyMatrix(s, refbasis)
    return CoherencyMatrix(t1*c_cir*t2, b1, b2)
end

function CoherencyMatrix(s::StokesParams, b1::T, b2::T, refbasis=CirBasis()) where {T<:PolBasis}
    return CoherencyMatrix(s, b1)
end

@inline function CoherencyMatrix(s::StokesParams, b::PolBasis)
    return CoherencyMatrix{typeof(b), typeof(b)}(s)
end

@inline function CoherencyMatrix{CirBasis,CirBasis}(s::StokesParams)
    (;I,Q,U,V) = s
    RR = complex((I + V))
    LR = (Q - 1im*U)
    RL = (Q + 1im*U)
    LL = complex((I - V))
    return CoherencyMatrix(RR, LR, RL, LL, CirBasis(), CirBasis())
end

@inline function CoherencyMatrix{B1, B2}(s::StokesParams) where {B1, B2}
    return CoherencyMatrix(s, B1(), B2())
end


@inline function CoherencyMatrix{LinBasis, LinBasis}(s::StokesParams)
    (;I,Q,U,V) = s
    XX = (I + Q)
    YX = (U - 1im*V)
    XY = (U + 1im*V)
    YY = (I - Q)
    return CoherencyMatrix(XX, YX, XY, YY, LinBasis(), LinBasis())
end



@inline function StokesParams(c::CoherencyMatrix{CirBasis, CirBasis})
    I = (c.e11 + c.e22)/2
    Q = (c.e21 + c.e12)/2
    U = 1im*(c.e21 - c.e12)/2
    V = (c.e11 - c.e22)/2
    return StokesParams(I, Q, U, V)
end

@inline function StokesParams(c::CoherencyMatrix{B1, B2}) where {B1, B2}
    t1 = basis_transform(B1()=>CirBasis())
    # Flip because these are the dual elements
    t2 = basis_transform(CirBasis()=>B2())
    c_cir = CoherencyMatrix(t1*c*t2, CirBasis())
    return StokesParams(c_cir)
end
