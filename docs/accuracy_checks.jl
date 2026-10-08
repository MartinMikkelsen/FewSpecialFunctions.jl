# Measurements behind the error figures quoted in docs/src/accuracy.md,
# docs/src/differentiation.md, and the Fermi–Dirac docstrings.
#
# Run from the repository root with
#
#     julia --project=docs docs/accuracy_checks.jl
#
# Every figure is a relative error |computed - reference| / |reference| unless
# labeled otherwise. References are closed forms, convergent series, or
# adaptive quadrature evaluated in 512-bit BigFloat arithmetic.

using FewSpecialFunctions, ForwardDiff, QuadGK, Printf

const REFERENCE_BITS = 512

relerr(computed, reference) = Float64(abs(big(computed) - reference) / abs(reference))

# Integral of `f` over `[lo, hi]` in high precision.
function reference_integral(f, lo, hi)
    return setprecision(BigFloat, REFERENCE_BITS) do
        first(quadgk(f, big(lo), big(hi); rtol = big"1e-40", order = 31))
    end
end

report(label, value) = @printf("  %-58s %s\n", label, value isa AbstractFloat ? @sprintf("%.2e", value) : string(value))

# ── Fermi–Dirac integrals ─────────────────────────────────────────────────────

# Unnormalized integral ∫₀^∞ t^j / (exp(t - x) + 1) dt. The substitution
# t = s² removes the endpoint singularity at j = -1/2; the integrand is below
# exp(-120) beyond the upper limit.
function fermi_dirac_reference(j, x)
    return setprecision(BigFloat, REFERENCE_BITS) do
        jb, xb = big(j), big(x)
        integrand(s) = 2 * s^(2jb + 1) / (exp(s^2 - xb) + 1)
        split = sqrt(max(xb, big(1)))
        upper = sqrt(max(xb, big(0)) + 120)
        first(quadgk(integrand, big(0), split, upper; rtol = big"1e-22", order = 21))
    end
end

const FERMI_GRID = vcat(-30.0:1.0:-3.0, -2.9:0.1:8.0, 9.0:1.0:60.0)
const FERMI_RATIONAL_ORDERS = (-1 // 2, 1 // 2, 3 // 2, 5 // 2)
const FERMI_OTHER_ORDERS = (1 // 4, 1 // 1, 2 // 1, 3 // 1, 9 // 2, 10 // 1)

function fermi_dirac_maximum(j)
    errors = [relerr(FermiDiracIntegral(Float64(j), x), fermi_dirac_reference(j, x)) for x in FERMI_GRID]
    i = argmax(errors)
    return errors[i], FERMI_GRID[i]
end

println(
    "Fermi–Dirac integrals: maximum relative error over $(length(FERMI_GRID)) points in ",
    "[$(first(FERMI_GRID)), $(last(FERMI_GRID))], Float64"
)
for j in (FERMI_RATIONAL_ORDERS..., 0 // 1, FERMI_OTHER_ORDERS...)
    err, x = fermi_dirac_maximum(j)
    report("order $(Float64(j))", @sprintf("%.2e at x = %g", err, x))
end
let grid = filter(>=(-5), FERMI_GRID)
    err = maximum(relerr(FermiDiracIntegral(0.0, x), fermi_dirac_reference(0, x)) for x in grid)
    report("order 0.0, restricted to x ≥ -5", err)
end
println("Fermi–Dirac integrals: 256-bit BigFloat arguments at x = 1")
for j in (3 // 2, 1 // 1)
    value = setprecision(() -> FermiDiracIntegral(big(j), big(1)), BigFloat, 256)
    report("order $(Float64(j))", relerr(value, fermi_dirac_reference(j, 1)))
end

# ── Spot checks of the other families ─────────────────────────────────────────

println("Coulomb wave functions at ℓ = 0, η = 0.3, ρ = 2")
let ℓ = 0, η = 3 // 10, ρ = 2
    at(f, T) = f(T(ℓ), T(η), T(ρ))
    high(f) = setprecision(() -> at(f, BigFloat), BigFloat, REFERENCE_BITS)
    report("F, Float64 vs 512-bit evaluation", relerr(at(F, Float64), high(F)))
    report("G, Float64 vs 512-bit evaluation", relerr(at(G, Float64), high(G)))
    report("F, 256-bit vs 512-bit evaluation", relerr(setprecision(() -> at(F, BigFloat), BigFloat, 256), high(F)))
end

println("Coulomb ℓ-derivatives at ℓ = 1, η = 0.3, ρ = 2 (default step)")
let η = 3 // 10, ρ = 2
    # Central difference with a step far below the default one, in 1024-bit arithmetic.
    function reference(f)
        return setprecision(BigFloat, 1024) do
            h = big"1e-120"
            (f(1 + h, big(η), big(ρ)) - f(1 - h, big(η), big(ρ))) / (2h)
        end
    end
    relabs(a, b) = Float64(abs(a - b) / abs(b))
    for (name, f) in (("Φ_dot", Φ), ("F_dot", F)), T in (Float64, BigFloat)
        dot = name == "Φ_dot" ? Φ_dot : F_dot
        value = setprecision(() -> dot(T(1), T(η), T(ρ)), BigFloat, 256)
        report("$name, $T order", relabs(big(value), reference(f)))
    end
    report("eps(Float64)^(2/3)", eps(Float64)^(2 / 3))
    report("eps(BigFloat)^(2/3) at 256 bits", Float64(setprecision(() -> eps(BigFloat)^(2 / big(3)), BigFloat, 256)))
end

println("Whittaker functions: W(0, 1/2, 2) = exp(-1), M(0, 1/2, 2) = 2 sinh(1)")
for T in (Float64, BigFloat)
    report("W, $T", relerr(WhittakerW(T(0), T(1) / 2, T(2)), exp(big(-1))))
    report("M, $T", relerr(WhittakerM(T(0), T(1) / 2, T(2)), 2 * sinh(big(1))))
end

println("Parabolic cylinder: D₁(2) = 2 exp(-1)")
for T in (Float32, Float64, BigFloat)
    report("ParabolicCylinderD, $T", relerr(ParabolicCylinderD(T(1), T(2)), 2 * exp(big(-1))))
end

println("Debye function at n = 3, β = 1, x = 2, default tolerance")
let reference = 3 / big(2)^3 * reference_integral(t -> t^3 / (exp(t) - 1), 0, 2)
    for T in (Float32, Float64, BigFloat)
        report("debye_function, $T", relerr(debye_function(T(3), T(1), T(2)), reference))
    end
end

println("Fresnel cosine integral at x = 1.5")
let reference = reference_integral(t -> cos(big(π) * t^2 / 2), 0, 3 // 2)
    for T in (Float32, Float64, BigFloat)
        report("FresnelC, $T", relerr(FresnelC(T(3) / 2), reference))
    end
end

println("Fresnel integrals: real BigFloat arguments that are evaluated")
for f in (FresnelC, FresnelS, FresnelE)
    status = map((3 // 2, 2, 3, 4, 5, 11 // 2, 6, 10, 40)) do x
        try
            f(BigFloat(x))
            "$(Float64(x)) ok"
        catch err
            err isa MethodError || rethrow()
            "$(Float64(x)) MethodError"
        end
    end
    report(string(f), join(status, ", "))
end

println("Dawson integral at x = 2.5")
let reference = exp(-big(5 // 2)^2) * reference_integral(t -> exp(t^2), 0, 5 // 2)
    for T in (Float32, Float64, BigFloat)
        report("dawson, $T", relerr(dawson(T(5) / 2), reference))
    end
end

println("Clausen function: Cl₂(π/2) = Catalan's constant")
let reference = setprecision(() -> big(Base.MathConstants.catalan), BigFloat, REFERENCE_BITS)
    report("Clausen(2, π/2), Float32", relerr(Clausen(2, Float32(π) / 2), reference))
    report("Clausen(2, π/2), Float64", relerr(Clausen(2, π / 2), reference))
    report("Clausen(2, π/2), BigFloat, N = 10", relerr(Clausen(2, big(π) / 2), reference))
    report("Clausen(2, π/2), BigFloat, N = 20", relerr(Clausen(2, big(π) / 2; N = 20), reference))
end

println("Bose–Einstein integral: B₁(-1) = Li₂(exp(-1))")
let reference = setprecision(() -> sum(exp(-big(n)) / big(n)^2 for n in 1:400), BigFloat, REFERENCE_BITS)
    for T in (Float32, Float64, BigFloat)
        report("BoseEinsteinIntegralNorm, $T", relerr(BoseEinsteinIntegralNorm(T(1), T(-1)), reference))
    end
end

println("Marcum Q-function at M = 1, a = 1.5")
let a = big(3 // 2)
    # Q₁(a, b) = ∫_b^∞ x exp(-(x² + a²)/2) I₀(a x) dx, with I₀ from its power series.
    besseli0(z) = sum((z / 2)^(2k) / factorial(big(k))^2 for k in 0:200)
    integrand(x) = x * exp(-(x^2 + a^2) / 2) * besseli0(a * x)
    for b in (3, 9)
        reference = reference_integral(integrand, b, 45)
        report("b = $b: Q = $(@sprintf("%.2e", reference)), Float32", relerr(MarcumQ(1.0f0, 1.5f0, Float32(b)), reference))
        report("b = $b, Float64", relerr(MarcumQ(1.0, 1.5, Float64(b)), reference))
        report("b = $b, BigFloat", relerr(MarcumQ(big(1), a, big(b)), reference))
    end
end

println("Voigt function at x = 0.7, y = 0.5")
let x = big(7 // 10), y = big(1 // 2)
    reference = reference_integral(t -> y / big(π) * exp(-t^2) / ((x - t)^2 + y^2), -14, 14)
    for T in (Float32, Float64, BigFloat)
        report("voigt, $T", relerr(voigt(T(7) / 10, T(1) / 2), reference))
    end
end

# ── Finite-difference derivatives of the ForwardDiff extension ────────────────

println("ForwardDiff finite-difference rules, relative error of the derivative")
let d = ForwardDiff.derivative
    # d𝓕_j/dx = j 𝓕_{j-1}, with the reference integral on the right-hand side.
    fermi_slope(x) = 3 // 2 * fermi_dirac_reference(1 // 2, x)
    report("FermiDiracIntegral(3/2, x) in x at x = 1", relerr(d(x -> FermiDiracIntegral(1.5, x), 1.0), fermi_slope(1)))
    report("FermiDiracIntegral(3/2, x) in x at x = 2", relerr(d(x -> FermiDiracIntegral(1.5, x), 2.0), fermi_slope(2)))
    report("FermiDiracIntegral(3/2, x) in x at x = 1, Float32", relerr(d(x -> FermiDiracIntegral(1.5f0, x), 1.0f0), fermi_slope(1)))
    # The Wronskian F'G - FG' equals 1.
    wronskian = d(ρ -> F(0.0, 0.3, ρ), 2.0) * G(0.0, 0.3, 2.0) - F(0.0, 0.3, 2.0) * d(ρ -> G(0.0, 0.3, ρ), 2.0)
    report("Coulomb Wronskian F'G - FG' - 1 (absolute)", abs(wronskian - 1))
    # dCl₁/dθ = -cot(θ/2)/2.
    report("Clausen(1, θ) in θ at θ = 1", relerr(d(θ -> Clausen(1, θ), 1.0), -cot(big(1) / 2) / 2))
    # ∂Q₁/∂a = b exp(-(a² + b²)/2) I₁(a b).
    besseli1(z) = sum((z / 2)^(2k + 1) / (factorial(big(k)) * factorial(big(k + 1))) for k in 0:200)
    a, b = big(3 // 2), big(3)
    report("MarcumQ(1, a, 3) in a at a = 1.5", relerr(d(a -> MarcumQ(1.0, a, 3.0), 1.5), b * exp(-(a^2 + b^2) / 2) * besseli1(a * b)))
end
