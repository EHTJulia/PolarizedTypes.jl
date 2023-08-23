"""
    basis_components([T=Float64,], e::ElectricFieldBasis, b::PolBasis)

Returns a static vector that contains the components of the electric field basis vector `e`
in terms of the polarization basis `b`. The first argument is optionally the eltype of the
static vector.

# Examples
```julia
julia> basis_components(Float64, R(), PolBasis{XPol,Y}())
2-element StaticArraysCore.SVector{2, ComplexF64} with indices SOneTo(2):
 0.7071067811865475 + 0.0im
                0.0 - 0.7071067811865475im

julia> basis_components(R(), PolBasis{XPol,Y}())
2-element StaticArraysCore.SVector{2, ComplexF64} with indices SOneTo(2):
 0.7071067811865475 + 0.0im
                0.0 - 0.7071067811865475im


julia> basis_components(Float64, X(), PolBasis{XPol,Y}())
2-element StaticArraysCore.SVector{2, ComplexF64} with indices SOneTo(2):
 1.0 + 0.0im
 0.0 + 0.0im
```
"""
function basis_components end


"""
    innerprod(::Type{T}, XPol(), YPol())

Computes the complex inner product of two elements of a complex Hilbert space `X` and `Y`
where base element of the output is T.
"""
function innerprod end

# Define that XPol,YPol and RPol,LPol are orthonormal bases
@inline innerprod(::Type{T}, ::B, ::B) where {T, B<:ElectricFieldBasis} = one(T)
@inline innerprod(::Type{T}, ::RPol, ::LPol) where {T} = complex(zero(T))
@inline innerprod(::Type{T}, ::LPol, ::RPol) where {T} = complex(zero(T))
@inline innerprod(::Type{T}, ::XPol, ::YPol) where {T} = complex(zero(T))
@inline innerprod(::Type{T}, ::YPol, ::XPol) where {T} = complex(zero(T))

# Now define the projections of linear onto circular
@inline innerprod(::Type{T}, ::XPol, ::RPol) where {T} = complex(inv(sqrt(T(2))), zero(T))
@inline innerprod(::Type{T}, ::XPol, ::LPol) where {T} = complex(inv(sqrt(T(2))), zero(T))

@inline innerprod(::Type{T}, ::YPol, ::RPol) where {T} = complex(zero(T), -inv(sqrt(T(2))))
@inline innerprod(::Type{T}, ::YPol, ::LPol) where {T} = complex(zero(T), inv(sqrt(T(2))))

# use the conjugate symmetry of the inner product to define projections of circular onto linear.
@inline innerprod(::Type{T}, c::Union{RPol,LPol}, l::Union{XPol,YPol}) where {T} = conj(innerprod(T, l, c))

# Now handle missing basis vectors (important when you are missing a feed)
@inline innerprod(::Type{T}, c::Missing, l::ElectricFieldBasis) where {T} = missing
@inline innerprod(::Type{T}, l::ElectricFieldBasis, c::Missing) where {T} = missing

# Now give me the components of electic fields in both linear and circular bases.
@inline basis_components(::Type{T}, b1::Union{ElectricFieldBasis,Missing}, ::PolBasis{B1,B2}) where {T, B1,B2} = SVector{2}(innerprod(T, B1(), b1), innerprod(T, B2(), b1))
@inline basis_components(v::Union{ElectricFieldBasis, Missing}, b::PolBasis) = basis_components(Float64, v, b)

#This handles non-orthogonal bases
# @inline basis_components(::Type{T}, b1::B1, ::PolBasis{B1, B2}) where {T, B1<:ElectricFieldBasis, B2<:ElectricFieldBasis} = SVector{2}(complex(one(T)), complex(zero(T)))
# @inline basis_components(::Type{T}, b1::B2, ::PolBasis{B1, B2}) where {T, B1<:ElectricFieldBasis, B2<:ElectricFieldBasis} = SVector{2}(complex(zero(T)), complex(one(T)))
# @inline basis_components(::Type{T}, b1::B1, ::PolBasis{B1, B1}) where {T, B1<:ElectricFieldBasis} = SVector{2}(complex(one(T)), complex(one(T)))


for (E1, E2) in [(:XPol,:RPol), (:RPol, :XPol), (:RPol,:YPol), (:YPol,:RPol), (:XPol,:LPol), (:LPol,:XPol), (:YPol,:LPol), (:LPol,:YPol)]
    @eval begin
        PolBasis{$E1,$E2}() = throw(AssertionError("Non-orthogonal bases not implemented"))
    end
end

"""
    basis_transform([T=Float64,], b1::PolBasis, b2::PolBasis)
    basis_transform([T=Float64,], b1::PolBasis=>b2::PolBasis)

Produces the transformation matrix that transforms the vector components from basis `b1` to basis `b2`.
This means that if for example `E` is the circular basis then `basis_transform(CirBasis=>LinBasis)E` is in the
linear basis. In other words the **columns** of the transformation matrix are the coordinate vectors
of the new basis vectors in the old basis.

# Example
```julia-repl
julia> basis_transform(CirBasis()=>LinBasis())
2×2 StaticArraysCore.SMatrix{2, 2, ComplexF64, 4} with indices SOneTo(2)×SOneTo(2):
 0.707107-0.0im       0.707107-0.0im
      0.0-0.707107im       0.0+0.707107im
```
"""
function basis_transform end

@inline basis_transform(::Type{T}, b1::PolBasis{E1,E2}, b2::PolBasis) where {E1,E2,T} = hcat(basis_components(T, E1(), b2), basis_components(T, E2(), b2))
# @inline basis_transform(::Type{T}, b1::B, b2::B) where {T, B<:PolBasis} = SMatrix{2,2,Complex{T}}(1.0, 0.0, 0.0, 1.0)
@inline basis_transform(b1::PolBasis, b2::PolBasis) = basis_transform(Float64, b1, b2)

@inline basis_transform(::Type{T}, p::Pair{B1,B2}) where {T, B1<:PolBasis, B2<:PolBasis} = basis_transform(T, B1(), B2())
@inline basis_transform(::Pair{B1,B2}) where {B1<:PolBasis, B2<:PolBasis} = basis_transform(Float64, B1(), B2())

@inline function basis_transform(c::CoherencyMatrix{B1, B2, Complex{T}}, e1::PolBasis, e2::PolBasis) where {B1, B2, T}
    t1 = basis_transform(T, B1()=>e1)
    t2 = basis_transform(T, e2=>B2())
    return CoherencyMatrix(t1*c*t2, e1, e2)
end
