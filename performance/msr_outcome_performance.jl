using LinearAlgebra
using SparseArrays
using QuantumOperators

function calc_1(m, v)
    return sum(abs2.(m * v))
end

function calc_2(m, v)
    out::Float64 = 0
    for j in axes(m, 2)
        out += abs2(dot(@view(m[:, j]), v)) 
    end
    return out
end

function f()
    d = 2; L = 18;
    msr_op = measurementoperators(nop(d), L)
    m = msr_op[1][1]
    v = rand(d^L)

    calc_1(m, v)
    calc_2(m, v)

    @time calc_1(m, v)
    @time calc_2(m, v)
end

f();