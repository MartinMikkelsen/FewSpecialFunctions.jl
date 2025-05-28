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