using LinearAlgebra

function calc_1(m, v)
    return sum(abs2.(m * v))
end

function calc_2(m, v)
    out::Float64 = 0
    for j in axes(m, 2)
        out += dot(v, @view(m[:, j])) 
    end
    return out::Float64
end

function f()
    n = 10000
    m = rand(n, n)
    v = rand(n)

    calc_1(m, v)
    calc_2(m, v)

    @time calc_1(m, v)
    @time calc_2(m, v)
end

f();