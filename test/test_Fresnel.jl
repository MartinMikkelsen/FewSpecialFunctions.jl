@testset "Fresnel" begin
    
    @testset "fresnel function" begin
        # Test that the fresnel function returns the correct values
        test_points = [0.0, 0.1, 0.5, 1.0, 2.0]
        
        for x in test_points
            # Test that the components match the individual functions
            S, C, E = fresnel(x)
            @test S ≈ FresnelS(x)
            @test C ≈ FresnelC(x)
            @test E ≈ FresnelE(x)
            
            # Test against reference implementations if available
            if isdefined(FewSpecialFunctions, :Fresnel_S_integral)
                @test S ≈ FewSpecialFunctions.Fresnel_S_integral(x) atol=1e-5
                @test C ≈ FewSpecialFunctions.Fresnel_C_integral(x) atol=1e-5
            end
        end
        
        # Test some known values
        S, C, _ = fresnel(0.0)
        @test S ≈ 0.0 atol=1e-10
        @test C ≈ 0.0 atol=1e-10
        
        # Test complex input
        S, C, E = fresnel(1.0 + 1.0im)
        @test isfinite(S) && isfinite(C) && isfinite(E)
        
        # Test that FresnelE behaves as expected
        for x in test_points
            S, C, E = fresnel(x)
            # E(z) = (1 + im)/2 - resp, where resp = (1 + im)/2 * erf((√π/2)*(1 + im)*z)
            # Could verify this relationship if needed
            @test isfinite(E)
        end
    end

end