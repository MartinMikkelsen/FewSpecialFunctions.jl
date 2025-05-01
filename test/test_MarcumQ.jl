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

    # Compared to marcumq (https://se.mathworks.com/help/signal/ref/marcumq.html)

    @test FewSpecialFunctions.marcum_Q(5, 0.47, 4.85) ≈ 0.0106766402997493 atol = 1e-9
    @test FewSpecialFunctions.marcum_Q(5, 1.46, 4.00) ≈ 0.211798804811782  atol = 1e-9
    @test FewSpecialFunctions.marcum_Q(1, 1.27, 4.58) ≈ 0.000931257801666407 atol = 1e-9
    @test FewSpecialFunctions.marcum_Q(4, 2.88, 3.28) ≈ 0.773155207859263 atol = 1e-9
    @test FewSpecialFunctions.marcum_Q(1, 2.55, 4.67) ≈  0.024112315799424 atol = 1e-9
    @test FewSpecialFunctions.marcum_Q(4, 2.27, 3.72) ≈ 0.400088995953665 atol = 1e-9
    @test FewSpecialFunctions.marcum_Q(2, 1.97, 0.86) ≈ 0.990345203236692 atol = 1e-9
    @test FewSpecialFunctions.marcum_Q(4, 0.10, 1.38) ≈ 0.983869651076909 atol = 1e-9
    @test FewSpecialFunctions.marcum_Q(1, 0.29, 4.12) ≈ 0.00028475422733874 atol = 1e-9
    @test FewSpecialFunctions.marcum_Q(4, 0.95, 4.75) ≈ 0.00899422673877906 atol = 1e-9

end