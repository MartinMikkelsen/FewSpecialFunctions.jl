@testset "Parabolic Cylinder function" begin
    
    # https://link.springer.com/content/pdf/10.1007/s00211-004-0517-x.pdf
    @test FewSpecialFunctions.U(10.1, 2*1.2*sqrt(10.1)) ≈ 8.7742145116891e-17 atol = 1e-9
    @test FewSpecialFunctions.U(20.1, 2*1.2*sqrt(20.1)) ≈ 2.8991030051243e-35 atol = 1e-9
    @test FewSpecialFunctions.U(30.1, 2*1.2*sqrt(30.1)) ≈  7.6172124886582e-55 atol = 1e-9

    @test FewSpecialFunctions.U(0.0, 10.0) ≈ exp(-0.25*10^2) atol = 1e-9
    @test FewSpecialFunctions.U(0.0, 20.0) ≈ exp(-0.25*20^2) atol = 1e-9
    @test FewSpecialFunctions.U(0.0, 30.0) ≈ exp(-0.25*30^2) atol = 1e-9

    @test FewSpecialFunctions.U(0.0, 40.0) ≈ exp(-0.25*40^2) atol = 1e-9
    @test FewSpecialFunctions.U(0.0, 50.0) ≈ exp(-0.25*50^2) atol = 1e-9
    @test FewSpecialFunctions.U(0.0, 60.0) ≈ exp(-0.25*60^2) atol = 1e-9
   

    # S. Zhang and J. Jin, 'Computation of Special functions' (Wiley, 1966),  E. Cojocaru, January 2009
    @test FewSpecialFunctions.U(-1.25459881152638, 5.70351922786027) ≈ 0.00109617508232108 atol = 1e-9
    @test FewSpecialFunctions.V(-1.25459881152638, 5.70351922786027) ≈ 139.15354241727 atol = 1e-9
    @test FewSpecialFunctions.W(-1.25459881152638, 5.70351922786027) ≈ 0.313946678917529 atol = 1e-9

    @test FewSpecialFunctions.U(4.50714306409916, -6.00652435683281) ≈ 1316297.50250584 atol = 1e-6
    @test FewSpecialFunctions.V(4.50714306409916, -6.00652435683281) ≈ 10162037.3095771 atol = 1e-6
    @test FewSpecialFunctions.W(4.50714306409916, -6.00652435683281) ≈ -28.5657116258641 atol = 1e-6

    @test FewSpecialFunctions.U(2.31993941811405, 0.284688768272233) ≈ 0.444692103073724 atol = 1e-9
    @test FewSpecialFunctions.V(2.31993941811405, 0.284688768272233) ≈ 0.784212080722314 atol = 1e-9
    @test FewSpecialFunctions.W(2.31993941811405, 0.284688768272233) ≈ 0.377302032929388 atol = 1e-9

    @test FewSpecialFunctions.dU(-1.25459881152638, 5.70351922786027) ≈ -0.00298204781259066 atol = 1e-9
    @test FewSpecialFunctions.dV(-1.25459881152638, 5.70351922786027) ≈ 349.325623311464 atol = 1e-9
    @test FewSpecialFunctions.dW(-1.25459881152638, 5.70351922786027) ≈ 1.41866642845807 atol = 1e-9

    @test FewSpecialFunctions.dU(4.50714306409916, -6.00652435683281) ≈  -4766982.76553318  atol = 1e-6
    @test FewSpecialFunctions.dV(4.50714306409916, -6.00652435683281) ≈ -36801905.8193446 atol = 1e-6
    @test FewSpecialFunctions.dW(4.50714306409916, -6.00652435683281) ≈ 2455.4339834424 atol = 1e-6

    @test FewSpecialFunctions.dU(2.31993941811405, 0.284688768272233) ≈  -0.693514373161243 atol = 1e-9
    @test FewSpecialFunctions.dV(2.31993941811405, 0.284688768272233) ≈ 0.57123166671827 atol = 1e-9
    @test FewSpecialFunctions.dW(2.31993941811405, 0.284688768272233) ≈-0.557336017048421 atol = 1e-9

end

