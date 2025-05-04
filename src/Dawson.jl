# --- Threshold accessor ---
function threshold(i)
    i == 1 && return 1.0e-12
    i == 2 && return 0.09
    i == 3 && return 0.3
    i == 4 && return 0.5
    i == 5 && return 0.7
    i == 6 && return 0.9
    i == 7 && return 1.1
    i == 8 && return 1.3
    i == 9 && return 30.0
    i == 10 && return 1.0e16
    error("Invalid threshold index")
end

# --- First 15 Chebyshev coefficients ---
const cffs = [
    2.46219053672763109e-35,
    2.46612209351226543e-33,
    2.53146774455005544e-31,
    2.66937551105311522e-29,
    2.89525328114429262e-27,
    3.23640440880034900e-25,
    3.73711810876911204e-23,
    4.46976231824321802e-21,
    5.55534876029332489e-19,
    7.20336188348558307e-17,
    9.79265074698045438e-15,
    1.40473769212007434e-12,
    2.14460623444015592e-10,
    3.52953377848361758e-08,
    6.38283985145170905e-06,
    1.33342567671828325e-03,
    1.46603629257834869e-03,
    1.91623518946637290e-16
]

# --- Continued fraction for small x ---
function priv_ncont_frac0(m, x)
    x2 = x^2
    coeffs = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]
    y = coeffs[m] * x2
    for k in (m - 1):-1:1
        y = coeffs[k] * x2 / (1 + y)
    end
    return x / (1 + y)
end

# --- Dawson integral ---
function daw_rk(x::BigFloat)
    ax = abs(x)

    if ax > threshold(10)
        return 0.5 / x
    elseif ax <= threshold(1)
        return x
    elseif ax <= threshold(2)
        return priv_ncont_frac0(1, x)
    elseif ax <= threshold(3)
        return priv_ncont_frac0(2, x)
    elseif ax <= threshold(4)
        return priv_ncont_frac0(3, x)
    elseif ax <= threshold(5)
        return priv_ncont_frac0(4, x)
    elseif ax <= threshold(6)
        return priv_ncont_frac0(5, x)
    elseif ax <= threshold(7)
        return priv_ncont_frac0(6, x)
    elseif ax <= threshold(8)
        return priv_ncont_frac0(7, x)
    else
        y100 = 100.0 / (ax + 1.8)
        ycase = floor(Int, y100)
        t = 2 * y100 - (2 * ycase + 1)
        res = t * cffs[1]
        for j in 2:14
            res = t * (res + cffs[j])
        end
        res = (res + cffs[15])
        return copysign(res, x)
    end
end

function daw_rk(x::Float64)
    return Float64(daw_rk(BigFloat(x)))
end
