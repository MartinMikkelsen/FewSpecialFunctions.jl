using SpecialFunctions

export FermiDiracIntegral

@doc raw"""
    FermiDiracIntegral(j,x)

The Fermi-Dirac integral
    
```math
    F_j(x) = \frac{1}{(\Gamma(j+1))} \int_0^\infty \frac{t^j}{\exp(t-x)+1} \, dt
```

Returns the value ``F_j(x)``

Resources:
[1] D. Bednarczyk and J. Bednarczyk, Phys. Lett. A, 64, 409 (1978)
[2] J. S. Blakemore, Solid-St. Electron, 25, 1067 (1982)
[3] X. Aymerich-Humet, F. Serra-Mestres, and J. Millan, Solid-St. Electron, 24, 981 (1981)
[4] X. Aymerich-Humet, F. Serra-Mestres, and J. Millan, J. Appl. Phys., 54, 2850 (1983)
https://arxiv.org/abs/0811.0116
https://de.wikipedia.org/wiki/Fermi-Dirac-Integral
https://dlmf.nist.gov/25.12#iii

"""
function FermiDiracIntegral(j::Real, eta::Real)
    if j < -1 / 2
        error("The order should be equal to or larger than -1/2.")
    else
        x = eta
        y = zero(x)
        if j == 0
            y = log(1 + exp(x))
        elseif j == 1 / 2
            # Model proposed in [1]
            # Expressions from eqs. (22)-(24) of [2]
            # mu = x^4 + 50 + 33.6 * x * (1 - 0.68 * exp(-0.17 * (x + 1)^2))
            # xi = 3 * sqrt(pi) / (4 * mu^(3 / 8))
            # y = (exp(-x) + xi)^-1
            a = 9.60
            b = 2.13
            c = 12/5
            y = (((j + 1) * 2^(j + 1) / (b + x + (abs(x - b).^ c + a) .^ (1 / c)) .^ (j + 1) +
                 exp(-x) / gamma(j + 1))^-1) / gamma(j + 1)
        elseif j == 3 / 2
            # Model proposed in [3]
            # Expressions from eq. (5) of [3]
            # The integral is divided by gamma(j + 1) to make it consistent with [1] and [2].
            a = 14.9
            b = 2.64
            c = 9 / 4
            y = (((j + 1) * 2^(j + 1) / (b + x + (abs(x - b).^ c + a) .^ (1 / c)) .^ (j + 1) +
                 exp(-x) / gamma(j + 1))^-1) / gamma(j + 1)
        else
            # Model proposed in [4]
            # Expressions from eqs. (6)-(7) of [4]
            # The integral is divided by gamma(j + 1) to make it consistent with [1] and [2].
            a = (1 + 15 / 4 * (j + 1) + 1 / 40 * (j + 1)^2)^(1 / 2)
            b = 1.8 + 0.61 * j
            c = 2 + (2 - sqrt(2)) * 2^(-j)
            y = ((j + 1) * 2^(j + 1) / (b + x + (abs(x - b) .^ c + a^c) .^ (1 / c)) .^ (j + 1) +
                 exp(-x) / gamma(j + 1))^-1 / gamma(j + 1)
        end
    end
    return y
end

