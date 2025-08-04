using PolarizedTypes
using StaticArrays
using ChainRulesCore
using JET
using Test

@testset "PolarizedTypes.jl" begin
    @testset "Basis transform" begin
        @test basis_transform(PolBasis{XPol,YPol}()=>PolBasis{RPol,LPol}())*basis_transform(PolBasis{RPol,LPol}()=>PolBasis{XPol,YPol}()) ≈ [1.0 0.0;0.0 1.0]
        @test basis_transform(PolBasis{RPol,LPol}()=>PolBasis{XPol,YPol}())*basis_transform(PolBasis{XPol,YPol}()=>PolBasis{RPol,LPol}()) ≈ [1.0 0.0;0.0 1.0]

        @test basis_transform(CirBasis(), LinBasis()) == basis_transform(CirBasis()=>LinBasis())

        @test basis_components(RPol(), CirBasis()) ≈ basis_transform(CirBasis()=>CirBasis())*SVector(1.0, 0.0)
        @test basis_components(LPol(), CirBasis()) ≈ basis_transform(CirBasis()=>CirBasis())*SVector(0.0, 1.0)
        @test basis_components(RPol(), LinBasis()) ≈ basis_transform(CirBasis()=>LinBasis())*SVector(1.0, 0.0)
        @test basis_components(LPol(), LinBasis()) ≈ basis_transform(CirBasis()=>LinBasis())*SVector(0.0, 1.0)

        @test basis_components(XPol(), CirBasis()) ≈ basis_transform(LinBasis()=>CirBasis())*SVector(1.0, 0.0)
        @test basis_components(YPol(), CirBasis()) ≈ basis_transform(LinBasis()=>CirBasis())*SVector(0.0, 1.0)
        @test basis_components(XPol(), LinBasis()) ≈ basis_transform(LinBasis()=>LinBasis())*SVector(1.0, 0.0)
        @test basis_components(YPol(), LinBasis()) ≈ basis_transform(LinBasis()=>LinBasis())*SVector(0.0, 1.0)


        for (e1, e2) in [(RPol, LPol), (LPol, RPol),
                         (XPol, YPol), (XPol, YPol),
                         ]
            @test basis_transform(PolBasis{e1,e2}()=>PolBasis{e1,e2}()) ≈ [1.0 0.0; 0.0 1.0]
        end
    end

    @testset "Non-orthogonal" begin
        @test_throws AssertionError basis_transform(PolBasis{XPol,YPol}(), PolBasis{RPol,XPol}())
        @test_throws AssertionError basis_transform(PolBasis{XPol,YPol}(), PolBasis{RPol,YPol}())
        @test_throws AssertionError basis_transform(PolBasis{XPol,YPol}(), PolBasis{XPol,RPol}())
        @test_throws AssertionError basis_transform(PolBasis{XPol,YPol}(), PolBasis{XPol,LPol}())
    end

    @testset "Missing feeds" begin
        for E in (XPol,YPol,RPol,LPol)
            c1 = basis_components(E(), PolBasis{E, Missing}())
            @test c1[1] ≈ 1.0
            @test c1[2] isa Missing

            c2 = basis_components(E(), PolBasis{Missing, E}())
            @test c2[2] ≈ 1.0
            @test c2[1] isa Missing

            # Test that the compiler realized the Union
            @inferred basis_components(E(), PolBasis{Missing, E}())
            @inferred basis_components(E(), PolBasis{E, Missing}())
            @inferred basis_transform(PolBasis{E, Missing}()=>CirBasis())
            @inferred basis_transform(PolBasis{E, Missing}()=>LinBasis())
            @inferred basis_transform(CirBasis()=>PolBasis{E, Missing}())
            @inferred basis_transform(LinBasis()=>PolBasis{E, Missing}())

        end
    end


    @testset "Simple stokes test" begin
        sQ = StokesParams(1.0, 0.5, 0.0, 0.0)
        sU = StokesParams(1.0, 0.0, 0.5, 0.0)
        sV = StokesParams(1.0, 0.0, 0.0, 0.5)

        @test CoherencyMatrix(sQ, LinBasis()) ≈ ([1.5 0.0; 0.0 0.5])
        @test CoherencyMatrix(sQ, CirBasis()) ≈ ([1.0 0.5; 0.5 1.0])

        @test CoherencyMatrix(sU, LinBasis()) ≈ ([1.0 0.5; 0.5 1.0])
        @test CoherencyMatrix(sU, CirBasis()) ≈ ([1.0 0.5im; -0.5im 1.0])

        @test CoherencyMatrix(sV, LinBasis()) ≈ ([1.0 0.5im; -0.5im 1.0])
        @test CoherencyMatrix(sV, CirBasis()) ≈ ([1.5 0.0; 0.0 0.5])
    end

    @testset "Simple Coherency test" begin
        cRR = CoherencyMatrix(0.5, 0.0, 0.0, 0.0, CirBasis())
        cLR = CoherencyMatrix(0.0, 0.5, 0.0, 0.0, CirBasis())
        cRL = CoherencyMatrix(0.0, 0.0, 0.5, 0.0, CirBasis())
        cLL = CoherencyMatrix(0.0, 0.0, 0.0, 0.5, CirBasis())

        @test StokesParams(cRR) ≈ inv(2)*[0.5, 0.0, 0.0, 0.5]
        @test StokesParams(cLR) ≈ inv(2)*[0.0, 0.5, 0.5im, 0.0]
        @test StokesParams(cRL) ≈ inv(2)*[0.0, 0.5, -0.5im, 0.0]
        @test StokesParams(cLL) ≈ inv(2)*[0.5, 0.0, 0.0, -0.5]


        cXX = CoherencyMatrix(0.5, 0.0, 0.0, 0.0, LinBasis())
        cYX = CoherencyMatrix(0.0, 0.5, 0.0, 0.0, LinBasis())
        cXY = CoherencyMatrix(0.0, 0.0, 0.5, 0.0, LinBasis())
        cYY = CoherencyMatrix(0.0, 0.0, 0.0, 0.5, LinBasis())

        @test StokesParams(cXX) ≈ inv(2)*[0.5, 0.5, 0.0, 0.0]
        @test StokesParams(cYX) ≈ inv(2)*[0.0, 0.0, 0.5, 0.5im]
        @test StokesParams(cXY) ≈ inv(2)*[0.0, 0.0, 0.5, -0.5im]
        @test StokesParams(cYY) ≈ inv(2)*[0.5, -0.5, 0.0, 0.0]

        @test StaticArraysCore.similar_type(CoherencyMatrix{CirBasis,LinBasis}, Float64, Size(2,2)) == CoherencyMatrix{CirBasis,LinBasis,Float64}
        @test StaticArraysCore.similar_type(StokesParams, Float64, Size(4,)) == StokesParams{Float64}


    end

    @testset "Conversions back and forward" begin
        s = StokesParams(1.0 .+ 0.0im, 0.2 + 0.2im, 0.2 - 0.2im, 0.1+0.05im)

        @test s ≈ StokesParams(CoherencyMatrix(s, CirBasis()))
        @test s ≈ StokesParams(CoherencyMatrix(s, LinBasis()))
        @test s ≈ StokesParams(CoherencyMatrix(s, CirBasis(), LinBasis()))
        @test s ≈ StokesParams(CoherencyMatrix(s, LinBasis(), CirBasis()))
        @test s ≈ StokesParams(CoherencyMatrix(s, PolBasis{YPol,XPol}(), PolBasis{LPol,RPol}()))
    end

    @testset "Mixed Pol" begin
        I = 2.0 + 0.5im
        Q = rand(ComplexF64) - 0.5
        U = rand(ComplexF64) - 0.5
        V = rand(ComplexF64) - 0.5
        s = StokesParams(I, Q, U, V)

        c1 = CoherencyMatrix(s, CirBasis(), LinBasis())
        c2 = CoherencyMatrix{CirBasis, LinBasis}(s)
        c3 = basis_transform(CoherencyMatrix(s, CirBasis()), CirBasis(), LinBasis())
        @test c1 ≈ c2
        @test c1 ≈ c3

        @test StokesParams(c1) ≈ s
        @test StokesParams(c2) ≈ s
        @test StokesParams(c3) ≈ s
    end

    @testset "Conversion Consistency" begin
        s = StokesParams(1.0 .+ 0.0im, 0.2 + 0.2im, 0.2 - 0.2im, 0.1+0.05im)
        c_lin1 = CoherencyMatrix(s, LinBasis())
        c_lin2 = CoherencyMatrix(s, PolBasis{XPol,YPol}())
        c_lin3 = CoherencyMatrix(s, PolBasis{XPol,YPol}(), PolBasis{XPol,YPol}())


        @test c_lin1 ≈ c_lin2 ≈ c_lin3

        c_cir1 = CoherencyMatrix(s, CirBasis())
        c_cir2 = CoherencyMatrix(s, PolBasis{RPol,LPol}())
        c_cir3 = CoherencyMatrix(s, PolBasis{RPol,LPol}(), PolBasis{RPol,LPol}())

        @test c_cir1 ≈ c_cir2 ≈ c_cir3

        t1 = basis_transform(LinBasis()=>CirBasis())
        t2 = basis_transform(CirBasis()=>LinBasis())

        @test t2*c_cir1*t1 ≈ c_lin1
        @test t1*c_lin1*t2 ≈ c_cir1

        @test_throws ArgumentError StokesParams(t1*c_lin1*t2)

        # Test the mixed basis
        @test c_cir1*t1 ≈ t1*c_lin1
        @test c_lin1*t2 ≈ t2*c_cir1

    end

    @testset "Performance test" begin
        s = StokesParams(1.0 .+ 0.0im, 0.2 + 0.2im, 0.2 - 0.2im, 0.1+0.05im)
        @test_opt StokesParams(CoherencyMatrix(s, LinBasis()))
        @test_opt StokesParams(CoherencyMatrix(s, CirBasis()))
        @test_opt StokesParams(CoherencyMatrix(s, LinBasis(), CirBasis()))
        @test_opt StokesParams(CoherencyMatrix(s, LinBasis(), CirBasis(), LinBasis()))
    end

    @testset "Polarized Functions" begin
        s = StokesParams(2.0, 0.25, -0.25, 0.25)
        @test linearpol(s) == complex(0.25, -0.25)
        @test evpa(s) ≈ atan(-0.5, 0.5)/2
        @test s[2:end] ≈ polarization(s)


        slin = StokesParams(2.0, 0.2, 0.2, 0.0)
        fp = fracpolarization(slin)
        @test complex(fp[1], fp[2]) ≈ linearpol(slin)/slin.I
        @test fp[end] ≈ 0

        @test fracpolarization(StokesParams(1f-5, 1f-5, 1f-6, 1f-7)) ≈ [1, 0.1, 0.01]

        @testset "ellipse" begin
            p = polellipse(s)
            @test p.a*p.b ≈ s.V^2/4
            @test p.evpa ≈ evpa(s)
            plin = polellipse(slin)
            @test plin.a ≈ (abs(linearpol(slin)))
            @test isapprox(plin.b, 0, atol=1e-8)
            @test plin.evpa ≈ evpa(slin)
            @test p.sn ≈ sign(s.V)
        end


        @test mpol(s) ≈ complex(0.25, -0.25)/2

        @testset "Complex Vis" begin
            I = 2.0 + 0.5im
            Q = rand(ComplexF64) - 0.5
            U = rand(ComplexF64) - 0.5
            V = rand(ComplexF64) - 0.5
            s = StokesParams(I, Q, U, V)
            @test mbreve(s) == m̆(s)
            c = CoherencyMatrix{CirBasis, LinBasis}(s)
            @test linearpol(c) ≈ linearpol(s)
            @test polarization(c) ≈ polarization(s)
            @test fracpolarization(c) ≈ fracpolarization(s)
            @test m̆(c) ≈ m̆(s)
            @test mbreve(c) ≈ mbreve(s)
            @test evpa(c) ≈ evpa(s)
        end

    end

    @testset "ChainRules" begin
        I = 2.0 + 0.5im
        Q = rand(ComplexF64) - 0.5
        U = rand(ComplexF64) - 0.5
        V = rand(ComplexF64) - 0.5
        s = StokesParams(I, Q, U, V)
        c = CoherencyMatrix{CirBasis, LinBasis}(s)

        cmat = SMatrix(c)
        prc =  ChainRulesCore.ProjectTo(c)
        @test prc(cmat) == c
        @test prc(c) == c
        @test_throws AssertionError prc(CoherencyMatrix(cmat, LinBasis()))
        @test_throws AssertionError prc(CoherencyMatrix(cmat, CirBasis()))
    end
end
