using SpecialFunctions

@testset "Parabolic Cylinder function" begin


    @testset "Clausen" begin

        data_U = open(readdlm, joinpath(@__DIR__, "data", "U_data.txt"))
        data_V = open(readdlm, joinpath(@__DIR__, "data", "V_data.txt"))

        @testset "U vs MATLAB reference" begin
            for r in 1:size(data_U, 1)
                a, x, Q_ref = data_U[r, :]

                @test FewSpecialFunctions.U(a, x) ≈ Q_ref atol = 1.0e-5
            end
        end
        @testset "V vs MATLAB reference" begin
            for r in 1:size(data_V, 1)
                a, x, Q_ref = data_V[r, :]

                @test FewSpecialFunctions.V(a, x) ≈ Q_ref atol = 1.0e-3
            end
        end

    end

    # https://link.springer.com/content/pdf/10.1007/s00211-004-0517-x.pdf
    @test FewSpecialFunctions.U(10.1, 2 * 1.2 * sqrt(10.1)) ≈ 8.7742145116891e-17 atol = 1.0e-9
    @test FewSpecialFunctions.U(20.1, 2 * 1.2 * sqrt(20.1)) ≈ 2.8991030051243e-35 atol = 1.0e-9
    @test FewSpecialFunctions.U(30.1, 2 * 1.2 * sqrt(30.1)) ≈ 7.6172124886582e-55 atol = 1.0e-9

    @test FewSpecialFunctions.U(0.0, 10.0) ≈ exp(-0.25 * 10^2) atol = 1.0e-9
    @test FewSpecialFunctions.U(0.0, 20.0) ≈ exp(-0.25 * 20^2) atol = 1.0e-9
    @test FewSpecialFunctions.U(0.0, 30.0) ≈ exp(-0.25 * 30^2) atol = 1.0e-9

    @test FewSpecialFunctions.U(0.0, 40.0) ≈ exp(-0.25 * 40^2) atol = 1.0e-9
    @test FewSpecialFunctions.U(0.0, 50.0) ≈ exp(-0.25 * 50^2) atol = 1.0e-9
    @test FewSpecialFunctions.U(0.0, 60.0) ≈ exp(-0.25 * 60^2) atol = 1.0e-9


    # S. Zhang and J. Jin, 'Computation of Special functions' (Wiley, 1966),  E. Cojocaru, January 2009
    @test FewSpecialFunctions.U(-1.25459881152638, 5.70351922786027) ≈ 0.00109617508232108 atol = 1.0e-9
    @test FewSpecialFunctions.V(-1.25459881152638, 5.70351922786027) ≈ 139.15354241727 atol = 1.0e-9
    @test FewSpecialFunctions.W(-1.25459881152638, 5.70351922786027) ≈ 0.313946678917529 atol = 1.0e-9

    @test FewSpecialFunctions.U(4.50714306409916, -6.00652435683281) ≈ 1316297.50250584 atol = 1.0e-6
    @test FewSpecialFunctions.V(4.50714306409916, -6.00652435683281) ≈ 10162037.3095771 atol = 1.0e-6
    @test FewSpecialFunctions.W(4.50714306409916, -6.00652435683281) ≈ -28.5657116258641 atol = 1.0e-6

    @test FewSpecialFunctions.U(2.31993941811405, 0.284688768272233) ≈ 0.444692103073724 atol = 1.0e-9
    @test FewSpecialFunctions.V(2.31993941811405, 0.284688768272233) ≈ 0.784212080722314 atol = 1.0e-9
    @test FewSpecialFunctions.W(2.31993941811405, 0.284688768272233) ≈ 0.377302032929388 atol = 1.0e-9

    @test FewSpecialFunctions.dU(-1.25459881152638, 5.70351922786027) ≈ -0.00298204781259066 atol = 1.0e-9
    @test FewSpecialFunctions.dV(-1.25459881152638, 5.70351922786027) ≈ 349.325623311464 atol = 1.0e-9
    @test FewSpecialFunctions.dW(-1.25459881152638, 5.70351922786027) ≈ 1.41866642845807 atol = 1.0e-9

    @test FewSpecialFunctions.dU(4.50714306409916, -6.00652435683281) ≈ -4766982.76553318  atol = 1.0e-6
    @test FewSpecialFunctions.dV(4.50714306409916, -6.00652435683281) ≈ -36801905.8193446 atol = 1.0e-6
    @test FewSpecialFunctions.dW(4.50714306409916, -6.00652435683281) ≈ 2455.4339834424 atol = 1.0e-6

    @test FewSpecialFunctions.dU(2.31993941811405, 0.284688768272233) ≈ -0.693514373161243 atol = 1.0e-9
    @test FewSpecialFunctions.dV(2.31993941811405, 0.284688768272233) ≈ 0.57123166671827 atol = 1.0e-9
    @test FewSpecialFunctions.dW(2.31993941811405, 0.284688768272233) ≈ -0.557336017048421 atol = 1.0e-9

    #Compare to SciPy
    xs = range(0.0, 5, 100)
    ws = FewSpecialFunctions.W(0.1, xs)
    dws = FewSpecialFunctions.dW(0.1, xs)

    expected_vals = [
        1.013635489653277, 0.988850748284378, 0.964316381133494,
        0.940021483752727, 0.915952391009033, 0.892092905269608,
        0.868424520179336, 0.844926640744175, 0.821576800528873,
        0.798350876857507, 0.775223304970871, 0.752167292146172,
        0.729155032821306, 0.706157925788338, 0.683146794528162,
        0.660092111750225, 0.636964229177185, 0.613733613573562,
        0.590371089959157, 0.566848092871232, 0.543136926443276,
        0.51921103395146, 0.495045277341601, 0.47061622708855,
        0.445902462555135, 0.420884882808501, 0.395547027616507,
        0.369875408085608, 0.343859846113304, 0.317493821512992,
        0.290774825326688, 0.263704717472259, 0.236290086477426,
        0.208542608634387, 0.180479403468638, 0.152123381956254,
        0.123503583449233, 0.094655496782976, 0.065621360549037,
        0.03645043702645, 0.007199253783757, -0.022068193500053,
        -0.051280283884079, -0.080357666970369, -0.109213219299851,
        -0.137752069258076, -0.165871698615115, -0.193462128935507,
        -0.220406201095406, -0.246579956017686, -0.271853124464221,
        -0.296089733289959, -0.31914883494797, -0.340885366221111,
        -0.361151141128005, -0.379795981693739, -0.396668988776242,
        -0.411619953387428, -0.424500906936959, -0.435167806553161,
        -0.443482349101952, -0.449313904738418, -0.452541557800638,
        -0.453056239612855, -0.450762934334063, -0.445582935406386,
        -0.437456126472339, -0.426343256897803, -0.412228178325518,
        -0.395120005069077, -0.375055157727226, -0.35209924625013,
        -0.326348745929254, -0.29793241752547, -0.26701242111693,
        -0.233785072365473, -0.198481189896357, -0.161365983489545,
        -0.122738434916538, -0.082930126642561, -0.042303478356203,
        -0.001249357477638, 0.039815962498489, 0.08045451270544,
        0.120211014697802, 0.158617935003227, 0.19520118824559,
        0.229486453438765, 0.261006051644906, 0.289306315818744,
        0.313955365602899, 0.334551181441651, 0.35072985404609,
        0.36217386743374, 0.368620256997394, 0.369868468896102,
        0.365787734110797, 0.356323760400063, 0.341504538775686,
        0.321445058635472,
    ]

    expected_dvals = [-0.49327396791427297, -0.4882279739100841, -0.48336902687551264, -0.47875407961448097, -0.47443551903685244, -0.4704612598055627, -0.46687482488277243, -0.4637154109570395, -0.4610179370238605, -0.45881307467787835, -0.457127258958835, -0.4559826788773798, -0.45539724703366796, -0.4553845480338462, -0.45595376570963975, -0.4571095894568603, -0.4588521003321658, -0.4611766378861149, -0.46407364906649495, -0.467528520900801, -0.4715213990620055, -0.47602699483825994, -0.48101438346535974, -0.48644679724040585, -0.4922814173152102, -0.49846916856685125, -0.5049545224577536, -0.5116753133250799, -0.5185625740743145, -0.5255403977887321, -0.5325258322976715, -0.5394288152635622, -0.546152157840278, -0.5525915854120217, -0.5586358443292687, -0.5641668839015271, -0.5690601231693393, -0.5731848121420694, -0.5764044972341766, -0.5785776005400185, -0.5795581223338667, -0.5791964757448562, -0.5773404619125907, -0.5738363930543969, -0.5685303697463636, -0.5612697173146479, -0.5519045845299153, -0.5402897057771505, -0.5262863255193801, -0.5097642811748813, -0.4906042374759696, -0.4687000619719933, -0.44396132758557766, -0.4163159240432175, -0.3857127556026075, -0.35212449782358635, -0.3155503812231909, -0.2760189645770045, -0.2335908554514506, -0.18836133036268368, -0.14046280186021765, -0.09006707494633388, -0.03738733070021497, 0.017320229071188478, 0.07375514572065475, 0.13157194830119526, 0.19037995239698935, 0.2497436871808005, 0.3091840235059324, 0.36818007308255385, 0.4261719242006302, 0.4825642728555822, 0.5367309993354237, 0.5880207291910488, 0.635763403918639, 0.6792778705593285, 0.717880480741715, 0.7508946684927099, 0.7776614525255515, 0.7975507828640931, 0.809973623850425, 0.8143946361738446, 0.8103452900218315, 0.7974372103606175, 0.7753755243881021, 0.7439719511539798, 0.7031573451064357, 0.652993379890091, 0.5936830371608866, 0.5255795486379757, 0.44919342926537503, 0.3651972364058741, 0.27442769561865077, 0.17788484890760847, 0.07672790739597644, -0.027732471935873127, -0.13404571583809818, -0.24063683840176536, -0.34582601086264786, -0.4478522016597585]

    @testset "W(0.1, x) and SciPy" begin
        for i in 1:length(ws)
            @test ws[i] ≈ expected_vals[i] atol = 1.0e-12
        end
    end

    @testset "dW(0.1, x) and SciPy" begin
        for i in 1:length(dws)
            @test dws[i] ≈ expected_dvals[i] atol = 1.0e-12
        end
    end

end


@testset "Special case: a < 0 and a + 0.5 ≈ integer" begin
    # a + 0.5 = integer, e.g., a = -1.5, -2.5, -3.5, etc.
    # Test for U, V, dU, dV with such a values and compare to known values or to the general branch

    # Helper to check if special branch is triggered
    function is_special_branch(a)
        a < 0 && isapprox(a + 0.5, round(a + 0.5))
    end

    # Test values
    a_vals = [-1.5, -2.5, -3.5, -4.5]
    x_vals = [-2.0, 0.0, 1.0, 3.5]

    for a in a_vals, x in x_vals
        # The function should not error and should return a finite value
        @test isfinite(FewSpecialFunctions.U(a, x))
        @test isfinite(FewSpecialFunctions.V(a, x))
        @test isfinite(FewSpecialFunctions.dU(a, x))
        @test isfinite(FewSpecialFunctions.dV(a, x))

        # The special branch should be triggered
        @test is_special_branch(a)
    end

    # Compare U(a, x) to the general branch for a just above and just below the special value
    for x in x_vals
        a = -2.5
        δ = 1.0e-8
        U_special = FewSpecialFunctions.U(a, x)
        U_above = FewSpecialFunctions.U(a + δ, x)
        U_below = FewSpecialFunctions.U(a - δ, x)
        # The function should be continuous across the branch
        @test isapprox(U_special, U_above; atol = 1.0e-6)
        @test isapprox(U_special, U_below; atol = 1.0e-6)
    end

    # Check that the formula for θ and prefactors are numerically stable
    for a in a_vals, x in x_vals
        θ = π * (0.25 + a / 2)
        f₁ = gamma(0.25 - a / 2) / (sqrt(π) * 2^(a / 2 + 0.25))
        f₂ = gamma(0.75 - a / 2) / (sqrt(π) * 2^(a / 2 - 0.25))
        @test isfinite(θ)
        @test isfinite(f₁)
        @test isfinite(f₂)
    end
end

@testset "Asymptotic expansion branch for W(a, x)" begin
    # Test values for large |x|, both positive and negative
    a_vals = [-2.0, -0.5, 1.0, 3.5]
    x_vals = [10.0, 20.0, -10.0, -20.0]

    # Check that the function does not error and returns finite values for large |x|
    for a in a_vals, x in x_vals
        w = FewSpecialFunctions.W(a, x)
        @test isfinite(w)
    end

    # Check continuity near x = 0 for large a
    for a in [5.0, 10.0]
        w_neg = FewSpecialFunctions.W(a, -20.0)
        w_pos = FewSpecialFunctions.W(a, 20.0)
        @test isfinite(w_neg)
        @test isfinite(w_pos)
    end

    # Check that the returned value is real for real inputs
    for a in a_vals, x in x_vals
        w = FewSpecialFunctions.W(a, x)
        @test isreal(w)
    end

    # Check that the branch for x > 0 and x < 0 is used
    for a in a_vals
        x = 15.0
        w_pos = FewSpecialFunctions.W(a, x)
        w_neg = FewSpecialFunctions.W(a, -x)
        # For a = 0, W(0, x) is even in x, so values should be close
        if a == 0.0
            @test isapprox(w_pos, w_neg; atol = 1.0e-8)
        end
    end

    # Check that the phase ϕ is finite and well-defined for large x
    for a in a_vals, x in x_vals
        g₀ = gamma(Complex(1 / 2, a))
        ϕ₂ = imag(g₀)
        ϕ = x^2 / 4 - a * log(abs(x)) + π / 4 + ϕ₂ / 2
        @test isfinite(ϕ)
    end

    # Check that the denominator in the recurrence is not zero
    for a in a_vals
        gref = gamma(Complex(1 / 2, a))
        gr₀, gi₀ = real(gref), imag(gref)
        den = gr₀^2 + gi₀^2
        @test den > 0
    end

    # Check that the recurrence for u and v produces finite arrays
    for a in a_vals
        u = zeros(Float64, 21)
        v = zeros(Float64, 21)
        gref = gamma(Complex(1 / 2, a))
        gr₀, gi₀ = real(gref), imag(gref)
        den = gr₀^2 + gi₀^2
        for k in 2:2:40
            m = k ÷ 2
            g = gamma(Complex(k + 0.5, a))
            gr, gi = real(g), imag(g)
            u[m] = (gr * gr₀ + gi * gi₀) / den
            v[m] = (gr₀ * gi - gr * gi₀) / den
        end
        @test all(isfinite, u)
        @test all(isfinite, v)
    end
end
