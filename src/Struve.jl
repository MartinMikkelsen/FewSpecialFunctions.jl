
using SpecialFunctions

export StruveH0, StruveH0Y0, StruveH1, StruveH1Y1

# Theo2 (2025). Struve functions (https://www.mathworks.com/matlabcentral/fileexchange/37302-struve-functions), MATLAB Central File Exchange. Retrieved July 21, 2025.

"""
    polyval(c::AbstractVector, x::AbstractVector)

Evaluate the polynomial with coefficients `c` at points `x` using Horner's method.
"""
function polyval(c::AbstractVector{<:Real}, x::AbstractVector{<:Real})
    result = fill(float(c[1]), length(x))
    for a in c[2:end]
        result .= result .* x .+ float(a)
    end
    return result
end

"""
    cheval(Ctype::String, An::AbstractVector, z::AbstractArray)

Evaluate Chebyshev series using the algorithm from Y.L. Luke (1977).
Supports the following types:
- `"regular"`: Tₙ(x)
- `"shifted"`: T*_n(x)
- `"even"`:    T_{2n}(x)
- `"odd"`:     T_{2n+1}(x)
"""
function cheval(Ctype::String, An::AbstractVector, z::AbstractArray)
    log1, log2 = 
        Ctype == "regular"  ? (true,  true)  :
        Ctype == "shifted"  ? (true,  false) :
        Ctype == "even"     ? (false, true)  :
        Ctype == "odd"      ? (false, false) :
        error("Unknown Chebyshev type: $Ctype")

    x = vec(z)
    xfac = 2 .* (2 .* x .- 1)
    if log1 && log2
        xfac = 2 .* x
    elseif !log1
        xfac = 2 .* (2 .* x.^2 .- 1)
    end

    n = length(An)
    Bn = zeros(eltype(x), length(x), n + 2)

    for j in 1:n
        Bn[:, end - j - 1] .= xfac .* Bn[:, end - j] .- Bn[:, end - j + 1] .+ An[n + 1 - j]
    end

    eval = Bn[:,1] .- xfac .* Bn[:,2] ./ 2
    if !log1 && !log2
        eval .= x .* (Bn[:,1] .- Bn[:,2])
    end

    return reshape(eval, size(z))
end

const nom_coeffs = Float64[
     4, 0, 8816, 0, 6778206,
     0, 2317961250, 0, 374638650750,
     0, 28306147182390, 0, 937687098467700,
     0, 11970864124436700, 0, 46174193179143750,
     0, 32071055725170000, 0, 840808761085125
]

const den_coeffs = Float64[
     4, 0, 8820, 0, 6786990,
     0, 2324669760, 0, 376904178000,
     0, 28663562736900, 0, 963414191990250,
     0, 12740661151218000, 0, 54025303535453250,
     0, 50023429199493750, 0, 4092826025413125
]

"""
    StruveH0Y0(z::AbstractArray)

Compute the function H₀(z) - Y₀(z) for complex arguments `z` using rational approximations for |z| > 16 and the exact form otherwise.
"""
function StruveH0Y0(z::AbstractArray)
    x = vec(z)
    f = similar(x, ComplexF64)

    i1 = abs.(x) .<= 16
    if any(i1)
        f[i1] .= StruveH0(x[i1]) .- bessely.(0, x[i1])
    end

    i2 = (abs.(x) .> 16) .& (real.(x) .< 0) .& (imag.(x) .< 0)
    if any(i2)
        x2 = -x[i2]
        P = polyval(nom_coeffs, x2)
        Q = polyval(den_coeffs, x2)
        f[i2] .= -(2 ./ π) ./ x2 .* (P ./ Q) .+ 2im .* besselh.(0, 1, x2)
    end

    i3 = (abs.(x) .> 16) .& (real.(x) .< 0) .& (imag.(x) .>= 0)
    if any(i3)
        x3 = -x[i3]
        P = polyval(nom_coeffs, x3)
        Q = polyval(den_coeffs, x3)
        f[i3] .= -(2 ./ π) ./ x3 .* (P ./ Q) .- 2im .* besselh.(0, 2, x3)
    end

    i4 = (abs.(x) .> 16) .& (real.(x) .>= 0)
    if any(i4)
        x4 = x[i4]
        P = polyval(nom_coeffs, x4)
        Q = polyval(den_coeffs, x4)
        f[i4] .= (2 ./ π) ./ x4 .* (P ./ Q)
    end

    return reshape(f, size(z))
end

const H0_bn = [
     7.741208505217204e-2,  -1.489049907224780e-1,   1.365908688119806e-1,
    -1.244408977167160e-1,   1.213210541746413e-1,  -9.902285991131098e-2,
     9.285469594924457e-2,  -8.889396057027860e-2,   6.318108596562945e-2,
    -3.164140591816664e-2,   1.160282350916835e-2,  -3.255171535152174e-3,
     7.240935462046915e-4,  -1.313122517655078e-4,   1.984234861137538e-5,
    -2.542464986238115e-6,   2.802187746192928e-7,  -2.688335553599889e-8,
     2.267665290294826e-9,  -1.696387437431367e-10,  1.133879020090459e-11,
    -6.816279772928244e-13,  3.706645666482425e-14, -1.832769429508806e-15,
     8.278232328724399e-17, -3.429926823742390e-18,  1.308568653462734e-19,
    -4.612908468868213e-21,  1.507287730295555e-22, -4.578534841587558e-24,
     1.296387199302720e-25, -3.430073819717729e-27,  8.500386329076394e-29,
    -1.977322978489543e-30,  4.326093533645215e-32
]

"""
    StruveH0(z::AbstractArray)

Evaluate the Struve function H₀(z) for complex `z` using a Chebyshev expansion for |z| ≤ 16 and asymptotics otherwise.
"""
function StruveH0(z::AbstractArray)
    x = vec(z)
    f = similar(x, ComplexF64)

    i1 = abs.(x) .<= 16
    if any(i1)
        z1 = x[i1].^2 ./ 400
        f[i1] .= cheval("shifted", H0_bn, z1) .* (2 .* x[i1]) ./ π
    end

    i2 = .!i1
    if any(i2)
        f[i2] .= StruveH0Y0(x[i2]) .+ bessely.(0, x[i2])
    end

    return reshape(f, size(z))
end

const H1_bn = [
    1.174772580755468e-1,  -2.063239340271849e-1,   1.751320915325495e-1,
   -1.476097803805857e-1,   1.182404335502399e-1,  -9.137328954211181e-2,
    6.802445516286525e-2,  -4.319280526221906e-2,   2.138865768076921e-2,
   -8.127801352215093e-3,   2.408890594971285e-3,  -5.700262395462067e-4,
    1.101362259325982e-4,  -1.771568288128481e-5,   2.411640097064378e-6,
   -2.817186005983407e-7,   2.857457024734533e-8,  -2.542050586813256e-9,
    2.000851282790685e-10, -1.404022573627935e-11,  8.842338744683481e-13,
   -5.027697609094073e-14,  2.594649322424009e-15, -1.221125551378858e-16,
    5.263554297072107e-18, -2.086067833557006e-19,  7.628743889512747e-21,
   -2.582665191720707e-22,  8.118488058768003e-24, -2.376158518887718e-25,
    6.492040011606459e-27, -1.659684657836811e-28,  3.978970933012760e-30,
   -8.964275720784261e-32,  1.901515474817625e-33
]

"""
    StruveH1(z::AbstractArray)

Evaluate the Struve function H₁(z) for complex arguments `z` using a Chebyshev expansion for |z| ≤ 16 and asymptotics otherwise.
"""
function StruveH1(z::AbstractArray)
    x = vec(z)
    f = similar(x, ComplexF64)

    i1 = abs.(x) .<= 16
    if any(i1)
        z1 = x[i1].^2 ./ 400
        f[i1] .= cheval("shifted", H1_bn, z1) .* x[i1].^2 .* (2 / (3π))
    end

    i2 = .!i1
    if any(i2)
        f[i2] .= StruveH1Y1(x[i2]) .+ bessely.(1, x[i2])
    end

    return reshape(f, size(z))
end


const H1Y1_nom = Float64[
     4, 0, 9648, 0, 8187030,
     0, 3.12092235e9, 0, 5.6883921003e11,
     0, 4.910820858405e13, 0, 1.8840528532161e15,
     0, 2.81319141807585e16, 0, 1.2623252631672375e17,
     0, 9.7007862050064e16, 0, 2.246438344775625e15
]

const H1Y1_den = Float64[
     4, 0, 9660, 0, 8215830,
     0, 3.14514144e9, 0, 5.779197396e11,
     0, 5.07124571499e13, 0, 2.01441149234325e15,
     0, 3.2559467386446e16, 0, 1.7751171161648925e17,
     0, 2.3010777431767125e17, 0, 3.1378332861500625e16
]

"""
    StruveH1Y1(z::AbstractArray)

Compute the function H₁(z) - Y₁(z) for complex arguments `z` using rational approximations for |z| > 16 and the exact form otherwise.
"""
function StruveH1Y1(z::AbstractArray)
    x = vec(z)
    f = similar(x, ComplexF64)

    i1 = abs.(x) .<= 16
    if any(i1)
        f[i1] .= StruveH1(x[i1]) .- bessely.(1, x[i1])
    end

    i2 = (abs.(x) .> 16) .& (real.(x) .< 0) .& (imag.(x) .< 0)
    if any(i2)
        x2 = -x[i2]
        num = polyval(H1Y1_nom, x2)
        den = polyval(H1Y1_den, x2)
        f[i2] .= 2/π .+ (2 ./ π) ./ x2.^2 .* (num ./ den) .- 2im .* besselh.(1, 1, x2)
    end

    i3 = (abs.(x) .> 16) .& (real.(x) .< 0) .& (imag.(x) .>= 0)
    if any(i3)
        x3 = -x[i3]
        num = polyval(H1Y1_nom, x3)
        den = polyval(H1Y1_den, x3)
        f[i3] .= 2/π .+ (2 ./ π) ./ x3.^2 .* (num ./ den) .+ 2im .* besselh.(1, 2, x3)
    end

    i4 = (abs.(x) .> 16) .& (real.(x) .>= 0)
    if any(i4)
        x4 = x[i4]
        num = polyval(H1Y1_nom, x4)
        den = polyval(H1Y1_den, x4)
        f[i4] .= 2/π .+ (2 ./ π) ./ x4.^2 .* (num ./ den)
    end

    return reshape(f, size(z))
end
