using LinearAlgebra
using SparseArrays
using QuantumOperators

function calc_1(m, v)
    return sum(abs2.(m * v))
end

function calc_2(m, v)
    out::Float64 = 0
    for j in axes(m, 2)
        for r in nzrange(m, j)
            i = rowvals(m)[r]
            out += abs2(nonzeros(m)[r]' * v[i]) 
        end
    end
    return out
end

function f()
    d = 2; L = 12;
    msr_op = measurementoperators(nop(d), L)
    m = msr_op[1][1]
    v = complex(rand(d^L))
    normalize!(v)

    calc_1(m, v)
    calc_2(m, v)

    @time calc_1(m, v)
    @time calc_2(m, v)

    measuresite!(v, msr_op, 1)
    @time measuresite!(v, msr_op, 1)
    # @time for i in 1:100 measuresite!(v, msr_op, 1) end
end

f();