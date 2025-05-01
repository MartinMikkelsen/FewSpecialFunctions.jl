using QuanticsTCI
using SpecialFunctions: besseli

export marcum_Q

function integrate_qtt(::Type{ValueType}, f::Function, a::Float64, b::Float64; N::Int=2^30, tolerance::Float64=1e-15) where {ValueType}
    xvals = range(a, b; length=N)
    dx = step(xvals)
    qtt, _, _ = quanticscrossinterpolate(ValueType, f, [xvals]; tolerance=tolerance)
    return dx * sum(qtt)
end

function marcum_Q(ν::Real, a::Real, b::Real)::Float64
    f(x) = (isa(x, Number) ? x : x[1])^ν * exp(-((isa(x, Number) ? x : x[1])^2 + a^2)/2) * besseli(ν - 1, a * (isa(x, Number) ? x : x[1]))
    I = 1 - 1 / (a^(ν - 1)) * integrate_qtt(Float64, f, 1e-18, b)
    return I
end

function marcum_Q(ν::Real, a::Real, b_array::AbstractArray{<:Real})::Vector{Float64}
    return [marcum_Q(ν, a, b_val) for b_val in b_array]
end

function marcum_Q(ν::Real, a_array::AbstractArray{<:Real}, b::Real)::Vector{Float64}
    return [marcum_Q(ν, a_val, b) for a_val in a_array]
end

function marcum_Q(ν_array::AbstractArray{<:Real}, a::Real, b::Real)::Vector{Float64}
    return [marcum_Q(ν_val, a, b) for ν_val in ν_array]
end
