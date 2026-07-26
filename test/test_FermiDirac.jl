@testset "Fermi-Dirac" begin

    # Comparing to https://npplus.readthedocs.io/en/latest/fermi.html
    # Test the output for j = -1/2
    @test FewSpecialFunctions.FermiDiracIntegral(-1 / 2, 0.0) ≈ 1.0721549299400754 atol = 1.0e-12
    @test FewSpecialFunctions.FermiDiracIntegral(-1 / 2, 1.0) ≈ 1.8204113571471041 atol = 1.0e-12
    @test FewSpecialFunctions.FermiDiracIntegral(-1 / 2, 1.2) ≈ 1.9785617633438695 atol = 1.0e-12
    @test FewSpecialFunctions.FermiDiracIntegral(-1 / 2, 4.1) ≈ 3.928454737099184 atol = 1.0e-12
    @test FewSpecialFunctions.FermiDiracIntegral(-1 / 2, 5.2) ≈ 4.477432715418454 atol = 1.0e-12
    @test FewSpecialFunctions.FermiDiracIntegral(-1 / 2, 6.6) ≈ 5.082787981164429 atol = 1.0e-12

    # Test the output for j = 1/2
    @test FewSpecialFunctions.FermiDiracIntegral(1 / 2, 0.0) ≈ 0.6780938951530457 atol = 1.0e-12
    @test FewSpecialFunctions.FermiDiracIntegral(1 / 2, 1.0) ≈ 1.3963752806666279 atol = 1.0e-12
    @test FewSpecialFunctions.FermiDiracIntegral(1 / 2, 1.2) ≈ 1.5863233997463857 atol = 1.0e-12
    @test FewSpecialFunctions.FermiDiracIntegral(1 / 2, 4.1) ≈ 5.965800008889902 atol = 1.0e-10
    @test FewSpecialFunctions.FermiDiracIntegral(1 / 2, 5.2) ≈ 8.28102917922544 atol = 1.0e-10
    @test FewSpecialFunctions.FermiDiracIntegral(1 / 2, 6.6) ≈ 11.632113406633252 atol = 1.0e-10

    # Test the output for j = 3/2
    @test FewSpecialFunctions.FermiDiracIntegral(3 / 2, 0.0) ≈ 1.15280383708879 atol = 1.0e-12
    @test FewSpecialFunctions.FermiDiracIntegral(3 / 2, 1.0) ≈ 2.6616826247307124 atol = 1.0e-12
    @test FewSpecialFunctions.FermiDiracIntegral(3 / 2, 1.2) ≈ 3.10869199517456 atol = 1.0e-12
    @test FewSpecialFunctions.FermiDiracIntegral(3 / 2, -2.5) ≈ 0.1075808743944384 atol = 1.0e-12
    @test FewSpecialFunctions.FermiDiracIntegral(3 / 2, -2.0) ≈ 0.175800988853926 atol = 1.0e-12

    # Test the output for j = 5/2
    @test FewSpecialFunctions.FermiDiracIntegral(5 / 2, 0.0) ≈ 3.0825860828379246 atol = 1.0e-13
    @test FewSpecialFunctions.FermiDiracIntegral(5 / 2, 1.0) ≈ 7.626535355004442 atol = 1.0e-13
    @test FewSpecialFunctions.FermiDiracIntegral(5 / 2, 1.2) ≈ 9.066754659807005 atol = 1.0e-13
    @test FewSpecialFunctions.FermiDiracIntegral(5 / 2, -2.5) ≈ 0.27085618652992494 atol = 1.0e-13
    @test FewSpecialFunctions.FermiDiracIntegral(5 / 2, -2.0) ≈ 0.4445544534586879 atol = 1.0e-13

    # Comparing to https://github.com/scott-maddox/fdint

    @test FewSpecialFunctions.FermiDiracIntegral(-1 / 2, 0.0) ≈ 1.0721549299401913 atol = 1.0e-12
    @test FewSpecialFunctions.FermiDiracIntegral(-1 / 2, 1.0) ≈ 1.8204113571469622 atol = 1.0e-12
    @test FewSpecialFunctions.FermiDiracIntegral(-1 / 2, 1.2) ≈ 1.9785617633438086 atol = 1.0e-12
    @test FewSpecialFunctions.FermiDiracIntegral(-1 / 2, 4.1) ≈ 3.928454737096438 atol = 1.0e-11
    @test FewSpecialFunctions.FermiDiracIntegral(-1 / 2, 5.2) ≈ 4.477432715417047 atol = 1.0e-11
    @test FewSpecialFunctions.FermiDiracIntegral(-1 / 2, 6.6) ≈ 5.082787981168395 atol = 1.0e-11

    # Test the output for j = 1/2
    @test FewSpecialFunctions.FermiDiracIntegral(1 / 2, 0.0) ≈ 0.6780938951531009 atol = 1.0e-12
    @test FewSpecialFunctions.FermiDiracIntegral(1 / 2, 1.0) ≈ 1.396375280666564 atol = 1.0e-12
    @test FewSpecialFunctions.FermiDiracIntegral(1 / 2, 1.2) ≈ 1.586323399746276 atol = 1.0e-12
    @test FewSpecialFunctions.FermiDiracIntegral(1 / 2, 4.1) ≈ 5.965800008886879 atol = 1.0e-10
    @test FewSpecialFunctions.FermiDiracIntegral(1 / 2, 5.2) ≈ 8.281029179226932 atol = 1.0e-10
    @test FewSpecialFunctions.FermiDiracIntegral(1 / 2, 6.6) ≈ 11.632113406638712 atol = 1.0e-10

    # Test the output for j = 3/2
    @test FewSpecialFunctions.FermiDiracIntegral(3 / 2, 0.0) ≈ 1.1528038370883613 atol = 1.0e-12
    @test FewSpecialFunctions.FermiDiracIntegral(3 / 2, 1.0) ≈ 2.6616826247320042 atol = 1.0e-10
    @test FewSpecialFunctions.FermiDiracIntegral(3 / 2, 1.2) ≈ 3.108691995173995 atol = 1.0e-12
    @test FewSpecialFunctions.FermiDiracIntegral(3 / 2, -2.5) ≈ 0.10758087439445904 atol = 1.0e-12
    @test FewSpecialFunctions.FermiDiracIntegral(3 / 2, -2.0) ≈ 0.1758009888540129 atol = 1.0e-12

    # Test the output for j = 5/2
    @test FewSpecialFunctions.FermiDiracIntegral(5 / 2, 0.0) ≈ 3.0825860828374183 atol = 1.0e-12
    @test FewSpecialFunctions.FermiDiracIntegral(5 / 2, 1.0) ≈ 7.626535355005596 atol = 1.0e-10
    @test FewSpecialFunctions.FermiDiracIntegral(5 / 2, 1.2) ≈ 9.066754659806907 atol = 1.0e-13
    @test FewSpecialFunctions.FermiDiracIntegral(5 / 2, -2.5) ≈ 0.2708561865299367 atol = 1.0e-13
    @test FewSpecialFunctions.FermiDiracIntegral(5 / 2, -2.0) ≈ 0.44455445345876304 atol = 1.0e-13

    # High-branch (x >= 2) reference values: -gamma(j+1)*polylog(j+1, -exp(x)),
    # mpmath at 50 digits. Guards the Antia [5] a2/b2 tables at design accuracy.
    @testset "FermiDirac high-branch reference values" begin
        refs = (
            (-1 / 2, (2.5953945832884784, 3.2852167828877117, 4.3832564347115706, 6.2971372445338478, 8.9349726661669702)),
            (1 / 2, (2.5024578260071403, 3.9769853540479774, 7.8379760572930966, 21.344471492355183, 59.812795370358027)),
            (3 / 2, (5.5372536750083451, 10.353714864761452, 27.80244621574838, 134.27015996313987, 726.56828396517523)),
            (5 / 2, (17.529419187384262, 36.932065413293268, 127.48954491322323, 1034.6842541815338, 10590.639176614387)),
        )
        for (j, expected) in refs
            for (x, ref) in zip((2.0, 3.0, 5.0, 10.0, 20.0), expected)
                @test FewSpecialFunctions.FermiDiracIntegral(j, x) ≈ ref rtol = 1.0e-10
            end
        end
    end

    @testset "FermiDirac invalid order" begin
        # Order j < -1/2 is not supported
        @test_throws ErrorException FewSpecialFunctions.FermiDiracIntegral(-1.0, 0.0)
        @test_throws ErrorException FewSpecialFunctions.FermiDiracIntegral(-0.6, 1.0)
    end
end
