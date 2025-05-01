@testset "Marcum Q-function" begin

    # Acta Univ. Sapientiae, Mathematica, 3, 1 (2011) 60–7
    @test FewSpecialFunctions.marcum_Q(1, 0.2, 0.6) ≈ 0.838249985438908 atol = 1e-9
    @test FewSpecialFunctions.marcum_Q(5, 0.2, 0.6) ≈ 0.999998670306184 atol = 1e-9
    @test FewSpecialFunctions.marcum_Q(7.7, 0.2, 0.6) ≈  0.999999999927717 atol = 1e-9

    @test FewSpecialFunctions.marcum_Q(1, 1.2, 1.6) ≈ 0.501536568390858 atol = 1e-9
    @test FewSpecialFunctions.marcum_Q(5, 1.2, 1.6) ≈ 0.994346394491553 atol = 1e-9
    @test FewSpecialFunctions.marcum_Q(7.7, 1.2, 1.6) ≈   0.999944937223540 atol = 1e-9

    @test FewSpecialFunctions.marcum_Q(1, 2.2, 2.6) ≈ 0.426794627821735 atol = 1e-9
    @test FewSpecialFunctions.marcum_Q(5, 2.2, 2.6) ≈ 0.929671935077756  atol = 1e-9
    @test FewSpecialFunctions.marcum_Q(7.7, 2.2, 2.6) ≈ 0.993735633182201 atol = 1e-9

    @test FewSpecialFunctions.marcum_Q(1, 2.2, 2.6) ≈ 0.426794627821735 atol = 1e-9
    @test FewSpecialFunctions.marcum_Q(5, 2.2, 2.6) ≈ 0.929671935077756  atol = 1e-9
    @test FewSpecialFunctions.marcum_Q(7.7, 2.2, 2.6) ≈ 0.993735633182201 atol = 1e-9


end