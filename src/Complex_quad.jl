using QuadGK

function complex_quadrature(func,a,b)
    function real_func(x)
        return real(func(x))
    end
    function imag_func(x)
        imag(func(x))
    end
    real_integral = quadgk(real_func,a,b)
    imag_integral = quadgk(imag_func,a,b)
    return real_integral[1] + 1im*imag_integral[1]
end

export complex_quadrature