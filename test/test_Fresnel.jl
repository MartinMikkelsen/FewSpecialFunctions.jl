using Test, DelimitedFiles, FewSpecialFunctions

@testset "Fresnel" begin

  @testset "Fresnel values" begin
    # read 5 columns: x, fsin, fcos, fexp_real, fexp_imag
    data = open(readdlm, joinpath(@__DIR__, "data", "FresnelF.txt"))

    for r in 1:size(data, 1)
      x          = data[r, 1]
      fsin_ref   = data[r, 2]
      fcos_ref   = data[r, 3]
      # Ensure complex parsing is robust
      fexp_str   = string(data[r, 4])
      fexp_ref   = parse(ComplexF64, replace(fexp_str, "i" => "im"))

      @test isapprox(FewSpecialFunctions.FresnelS(x), fsin_ref;   rtol=1e-2)
      @test isapprox(FewSpecialFunctions.FresnelC(x), fcos_ref;   rtol=1e-2)
    end
  end

end
@testset "fresnel function" begin
    # Test fresnel function directly
    x_values = [0.0, 1.0, 2.0, -1.5]
    
    for x in x_values
        C, S, E = FewSpecialFunctions.fresnel(x)
        
        # Test individual components match wrapper functions
        @test C ≈ FewSpecialFunctions.FresnelC(x)
        @test S ≈ FewSpecialFunctions.FresnelS(x)
        @test E ≈ FewSpecialFunctions.FresnelE(x)
        
        # Test that E = C + i*S relationship holds
        @test E ≈ C + im*S
    end
    
    # Test special case z=0
    C, S, E = FewSpecialFunctions.fresnel(0.0)
    @test C ≈ 0.0
    @test S ≈ 0.0
    @test E ≈ 0.0
    
    # Test complex input
    z = 1.0 + 1.0im
    C, S, E = FewSpecialFunctions.fresnel(z)
    @test E ≈ C + im*S
    
    # Test against known values
    # At z = 1, C(1) ≈ 0.7798934, S(1) ≈ 0.4382591
    C, S, E = FewSpecialFunctions.fresnel(1.0)
    @test isapprox(C, 0.7798934, rtol=1e-6)
    @test isapprox(S, 0.4382591, rtol=1e-6)
end

@testset "Fresnel edge and special cases" begin

  # Test at purely imaginary input
  z = 2.0im
  C, S, E = FewSpecialFunctions.fresnel(z)
  # E should equal C + im*S
  @test E ≈ C + im*S

  # Test at a small value (Taylor expansion regime)
  z = 1e-8
  C, S, E = FewSpecialFunctions.fresnel(z)
  @test isapprox(C, z, atol=1e-8)
  @test isapprox(S, (π/6)*z^3, atol=1e-24)

  # Test at a negative real value
  z = -2.0
  C, S, E = FewSpecialFunctions.fresnel(z)
  # Fresnel integrals are odd/even functions:
  @test isapprox(C, -FewSpecialFunctions.FresnelC(2.0), rtol=1e-6)
  @test isapprox(S, -FewSpecialFunctions.FresnelS(2.0), rtol=1e-6)

end
