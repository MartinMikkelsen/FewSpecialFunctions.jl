@testset "hypergeometric" begin
    
        # Comparing to https://mpmath.org/doc/current/functions/hypergeometric.html
        @test FewSpecialFunctions.confluent_hypergeometric_1F1(2.0,-1/3,3.25) ≈ -2815.956856924817275640248 atol = 1e-12
        @test FewSpecialFunctions.confluent_hypergeometric_1F1(2.0,-1/3,-3.25) ≈ -1.145036502407444445553107 atol = 1e-13
        @test FewSpecialFunctions.confluent_hypergeometric_1F1(1.0,3.0,1.0) ≈ 1.436563656918090470720575 atol = 1e-13
        @test FewSpecialFunctions.confluent_hypergeometric_1F1(1.0,3.0,2.0) ≈ 2.194528049465325113615214 atol = 1e-13
        @test FewSpecialFunctions.confluent_hypergeometric_1F1(1.0,3.0,3.0) ≈ 3.574563760708370609095229 atol = 1e-12
        @test FewSpecialFunctions.confluent_hypergeometric_1F1(1.0,3.0,4.0) ≈ 6.199768754143029884763783 atol = 1e-13
        @test FewSpecialFunctions.confluent_hypergeometric_1F1(1.0,3.0,5.0) ≈ 11.39305272820612827368925 atol = 1e-13
        @test FewSpecialFunctions.confluent_hypergeometric_1F1(1.0,3.0,6.0) ≈ 22.02382186070750681157707 atol = 1e-13
        @test FewSpecialFunctions.confluent_hypergeometric_1F1(1.0,3.0,7.0) ≈ 44.43400646646769792913144 atol = 1e-13
        @test FewSpecialFunctions.confluent_hypergeometric_1F1(1.0,3.0,8.0) ≈ 92.87368709505400858573725 atol = 1e-12
        @test FewSpecialFunctions.confluent_hypergeometric_1F1(1.0,3.0,9.0) ≈ 199.8292327796391113014814 atol = 1e-12
    
        @test FewSpecialFunctions.confluent_hypergeometric_1F1(3.0,5.0,1.0) ≈ 1.858145824525628472970349 atol = 1e-13
        @test FewSpecialFunctions.confluent_hypergeometric_1F1(3.0,5.0,2.0) ≈ 3.583584148395975340845641 atol = 1e-13
        @test FewSpecialFunctions.confluent_hypergeometric_1F1(3.0,5.0,3.0) ≈ 7.149127521416741218190458 atol = 1e-12
        @test FewSpecialFunctions.confluent_hypergeometric_1F1(3.0,5.0,4.0) ≈ 14.69947969682181724071851 atol = 1e-13
        @test FewSpecialFunctions.confluent_hypergeometric_1F1(3.0,5.0,5.0) ≈ 31.03765920246417864253961 atol = 1e-13
        @test FewSpecialFunctions.confluent_hypergeometric_1F1(3.0,5.0,6.0) ≈ 67.0714655821225204347312 atol = 1e-13
        @test FewSpecialFunctions.confluent_hypergeometric_1F1(3.0,5.0,7.0) ≈ 147.8838581136278992759039 atol = 1e-12
        @test FewSpecialFunctions.confluent_hypergeometric_1F1(3.0,5.0,8.0) ≈ 331.800010276129905586689 atol = 1e-12
        @test FewSpecialFunctions.confluent_hypergeometric_1F1(3.0,5.0,9.0) ≈ 755.7993238341921982500408 atol = 1e-12
    
        @test FewSpecialFunctions.confluent_hypergeometric_1F1(4.0,6.0,1.0) ≈ 1.977997739010048220736756 atol = 1e-13
        @test FewSpecialFunctions.confluent_hypergeometric_1F1(4.0,6.0,2.0) ≈ 4.027359752673374431923931 atol = 1e-13
        @test FewSpecialFunctions.confluent_hypergeometric_1F1(4.0,6.0,3.0) ≈ 8.41618195634263401010581 atol = 1e-12
        @test FewSpecialFunctions.confluent_hypergeometric_1F1(4.0,6.0,4.0) ≈ 17.99942188535757471190946 atol = 1e-13
        @test FewSpecialFunctions.confluent_hypergeometric_1F1(4.0,6.0,5.0) ≈ 39.28921294851610073770073 atol = 1e-13
        @test FewSpecialFunctions.confluent_hypergeometric_1F1(4.0,6.0,6.0) ≈ 87.31486279164030426724414 atol = 1e-13
        @test FewSpecialFunctions.confluent_hypergeometric_1F1(4.0,6.0,7.0) ≈ 197.129299604566250358567 atol = 1e-12
        @test FewSpecialFunctions.confluent_hypergeometric_1F1(4.0,6.0,8.0) ≈ 451.2631718666678540871648 atol = 1e-12
        @test FewSpecialFunctions.confluent_hypergeometric_1F1(4.0,6.0,9.0) ≈ 1045.691875021159141369882 atol = 1e-12
end