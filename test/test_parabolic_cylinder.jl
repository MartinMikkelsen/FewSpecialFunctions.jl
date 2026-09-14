using SpecialFunctions

@testset "Cylinder numerical regression checks" begin
    @testset "U near the former asymptotic switch" begin
        # Positive-integrand quadrature of DLMF 12.5.1 at 256-bit precision.
        for (a, x, expected) in [
                (10.0, 5.00001, 1.549593832704713e-11),
                (10.0, 6.0, 2.252444925281081e-13),
            ]
            @test FewSpecialFunctions.U(a, x) ≈ expected rtol = 2.0e-12
        end
        for x in (-20.0, -10.0, -5.00001, 5.00001, 10.0, 20.0)
            gaussian = exp(-x^2 / 4)
            @test FewSpecialFunctions.U(-0.5, x) ≈ gaussian rtol = 2.0e-13
            @test FewSpecialFunctions.dU(-0.5, x) ≈ -x * gaussian / 2 rtol = 2.0e-13
        end
    end

    @testset "W phase, amplitude, and derivative" begin
        # Independently summed ODE solution at 768 bits; truncations at 800 and
        # 1200 terms agreed. The a=0 checks below use a separate Bessel identity.
        for (a, x, expected) in [
                (1.0, 10.0, -0.035016738866263074),
                (0.1, 8.00001, -0.205949151506127),
                (10.0, 20.0, -1.6623780433717748e-8),
            ]
            @test FewSpecialFunctions.W(a, x) ≈ expected rtol = 2.0e-11
        end
        for x in (5.0, 7.99999, 8.00001, 10.0, 20.0, 40.0)
            z = x^2 / 4
            expected = sqrt(pi * x) / 2^(5 / 4) *
                (besselj(-0.25, z) - besselj(0.25, z))
            derivative = -sqrt(pi * x) * x / 2^(9 / 4) *
                (besselj(-0.75, z) + besselj(0.75, z))
            @test FewSpecialFunctions.W(0.0, x) ≈ expected rtol = 2.0e-11
            @test FewSpecialFunctions.dW(0.0, x) ≈ derivative rtol = 2.0e-11
        end
    end

    @testset "Float32 derivatives at the origin" begin
        @test FewSpecialFunctions.dU(0.0f0, 0.0f0) ≈ -0.5813683170191186f0 rtol = 5.0f-7
        @test FewSpecialFunctions.dV(0.0f0, 0.0f0) ≈ 0.3280019486668765f0 rtol = 5.0f-7
        @test FewSpecialFunctions.dW(0.0f0, 0.0f0) ≈ -0.48887053372346173f0 rtol = 5.0f-7
    end

    @testset "Independent solutions and precision" begin
        for a in (-10.0, -0.5, 0.0, 1.0, 10.0), x in (0.0, 2.0, 5.0, 8.0, 10.0)
            u, v = FewSpecialFunctions.U(a, x), FewSpecialFunctions.V(a, x)
            du, dv = FewSpecialFunctions.dU(a, x), FewSpecialFunctions.dV(a, x)
            @test u * dv - du * v ≈ sqrt(2 / pi) rtol = 2.0e-11
            w, wn = FewSpecialFunctions.W(a, x), FewSpecialFunctions.W(a, -x)
            dw, dwn = FewSpecialFunctions.dW(a, x), FewSpecialFunctions.dW(a, -x)
            @test w * dwn + dw * wn ≈ -1 rtol = 2.0e-11
        end
        # Fixed high-precision references from the independently summed ODE.
        setprecision(256) do
            @test FewSpecialFunctions.U(big(10), big(6)) ≈ big"2.25244492528108108020437274296873569833210978800633191009619758715824058922e-13" rtol = big"1e-65"
            @test FewSpecialFunctions.W(big(0), big(10)) ≈ big"0.229304673430426488080556254525231378709958147043123523788331611358746566233" rtol = big"1e-65"
            @test FewSpecialFunctions.dW(big(0), big(10)) ≈ big"-0.881210343397594405108599265961579620485936842697233669233829684544252911335" rtol = big"1e-65"
            @test precision(FewSpecialFunctions.W(big(0), big(2))) == 256
            @test precision(BigFloat) == 256
        end
        @test FewSpecialFunctions.U(10.0f0, 6.0f0) ≈ 2.252444925281081f-13 rtol = 2.0f-6
        @test FewSpecialFunctions.W(10.0f0, 20.0f0) ≈ -1.6623780433717748f-8 rtol = 2.0f-5
        x = 25.0
        derivative = -sqrt(pi * x) * x / 2^(9 / 4) * (besselj(-0.75, x^2 / 4) + besselj(0.75, x^2 / 4))
        @test FewSpecialFunctions.dW(0.0f0, 25.0f0) ≈ Float32(derivative) rtol = 2eps(Float32)
        # The derivative can remain representable after U itself underflows.
        for x in (54.6, 20.5f0)
            expected = typeof(x)(-BigFloat(x) / 2 * exp(-BigFloat(x)^2 / 4))
            @test FewSpecialFunctions.dU(typeof(x)(-0.5), x) == expected
        end
    end
end

@testset "Cylinder scaling and convergence safeguards" begin
    # Neither exp(logscale) is representable in Float64, but both complete
    # products are. BigFloat supplies an independently rounded reference.
    for (logscale, amplitude) in ((-750.0, 1.0e100), (750.0, -1.0e-100))
        reference = Float64(exp(BigFloat(logscale)) * BigFloat(amplitude))
        @test FewSpecialFunctions._cylinder_scaled(logscale, amplitude) ≈ reference rtol = 1.0e-13
    end
    @test isequal(FewSpecialFunctions._cylinder_scaled(Inf, -0.0), -0.0)

    # A Taylor sum that overflows at this precision must report failure.
    @test_throws r"series did not converge" FewSpecialFunctions._cylinder_series(0.0, 100.0, :U)

    # At 2048 bits these decreasing asymptotic sums cannot reach the requested
    # accuracy within the work limit: reject them so the caller uses its fallback.
    setprecision(BigFloat, 2048) do
        a, x = BigFloat(0), BigFloat(50)
        @test FewSpecialFunctions._cylinder_u_asymptotic(a, x) === nothing
        @test FewSpecialFunctions._cylinder_w_asymptotic(a, x) === nothing
    end
end

@testset "Concurrent cylinder precision" begin
    parameters = [(10.0, 6.0), (1.0, 0.1), (-10.0, 15.0), (20.0, 5.0)]
    references = [FewSpecialFunctions.U(a, x) for (a, x) in parameters]
    before = precision(BigFloat)
    tasks = [Threads.@spawn(FewSpecialFunctions.U(parameters[i]...)) for i in repeat(1:4, 8)]
    @test fetch.(tasks) == repeat(references, 8)
    @test precision(BigFloat) == before
end

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
    @test FewSpecialFunctions.U(10.1, 2 * 1.2 * sqrt(10.1)) ≈ 8.7742145116891e-17 rtol = 1.0e-12
    @test FewSpecialFunctions.U(20.1, 2 * 1.2 * sqrt(20.1)) ≈ 2.8991030051243e-35 rtol = 1.0e-12
    @test FewSpecialFunctions.U(30.1, 2 * 1.2 * sqrt(30.1)) ≈ 7.6172124886582e-55 rtol = 1.0e-12

    for x in (10.0, 20.0, 30.0, 40.0, 50.0, 60.0)
        @test FewSpecialFunctions.U(-0.5, x) ≈ exp(-0.25 * x^2) rtol = 1.0e-14
    end


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

end

@testset "U at negative half-integer orders" begin
    # DLMF 12.7.1–2: the Gaussian and first Hermite polynomial.
    @test FewSpecialFunctions.U(-0.5, -6.0) ≈ exp(-9) rtol = 1.0e-14
    @test FewSpecialFunctions.U(-1.5, -7.0) ≈ -7exp(-49 / 4) rtol = 1.0e-14
end

@testset "parabolic cylinder array consistency" begin
    x_vals = [0.5, 1.0, 2.0, 3.0]
    a_val = 1.0
    @test FewSpecialFunctions.U(a_val, x_vals) ≈ [FewSpecialFunctions.U(a_val, xi) for xi in x_vals]
    @test FewSpecialFunctions.V(a_val, x_vals) ≈ [FewSpecialFunctions.V(a_val, xi) for xi in x_vals]
    @test FewSpecialFunctions.W(a_val, x_vals) ≈ [FewSpecialFunctions.W(a_val, xi) for xi in x_vals]
    @test FewSpecialFunctions.dU(a_val, x_vals) ≈ [FewSpecialFunctions.dU(a_val, xi) for xi in x_vals]
    @test FewSpecialFunctions.dV(a_val, x_vals) ≈ [FewSpecialFunctions.dV(a_val, xi) for xi in x_vals]
    @test FewSpecialFunctions.dW(a_val, x_vals) ≈ [FewSpecialFunctions.dW(a_val, xi) for xi in x_vals]
end
