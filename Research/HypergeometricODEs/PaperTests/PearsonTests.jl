include("../pFq.jl")

tests = [
    (0.1,0.2,0.3,0.5)                   # 1
    (-0.1,0.2,0.3,0.5)                  # 2
    (.1,.2,-.3,-.5+.5im)                # 3
    (1e-8,1e-8,1e-8,1e-6)               # 4
    (1e-8,-1e-6,1e-12,-1e-10+1e-12im)   # 5

    (1,10,1,0.5+1e-9im)                 # 6
    (1,-1+1e-12im,1,-.8)                # 7
    (2+8im,3-5im,sqrt(2)-π*im,.75)      # 8
    (100,200,350,im)                    # 9
    (2+1e-9,3,5,-.75)                   # 10

    (-2,-3,-5+1e-9,.5)                  # 11
    (-1,-1.5,-2-1e-15,.5)               # 12
    (500,-500,500,.75)                  # 13
    (500,500,500,.75)                   # 14
    (-1000,-2000,-4000.1,-.5)           # 15

    (-100,-200,-300+1e-9,0.5sqrt(2))    # 16
    (300,10,5,.5)                       # 17
    (5,-300,10,.5)                      # 18
    (10,5,-300.5,.5)                    # 19
    (2+200im,5,10,0.6)                  # 20

    (2+200im,5-100im,10+500im,.8)       # 21
    (2,5,10-500im,-.8)                  # 22
    (2.25,3.75,-.5,-1)                  # 23
    (1,2,4+3im,.6-.8im)                 # 24
    (1,.9,2,cispi(1/3))                 # 25

    (1,1,4,cispi(1/3))                  # 26
    (-1,.9,2,cispi(-1/3))               # 27
    (4,1.1,2,.5+(.5sqrt(3)-0.01)*im)    # 28
    (5,2.2,-2.5,0.49+.5sqrt(3)*im)      # 29
    (2/3,1,4/3,cispi(1/3))              # 30
]

"Quick relative error test of 2f1."
function test_2f1(test::T; prec = 64) where T <: NTuple{4, Number}
    tru = mathematica_2f1(test...)

    if prec != 64
        setprecision(prec)
        val = taylor_2f1((big.(test))...)
    else
        val = taylor_2f1(test...)
    end

    return abs((val - tru) / tru)
end

"""
    run_pearson_tests(precs; order = 1000)
Run 2F1 tests mentioned in Pearson et al using the Taylor method of `order` and rerun test with extended `precs`ision calculations.
"""
function run_pearson_tests(precs = [128, 256, 512, 1024, 2048]; order = 5000)
    tru = [mathematica_2f1(test...) for test ∈ tests]

    vals = Matrix{Number}(undef, length(tests) + 1, length(precs) + 2)

    vals[1,:] = [[0, 64]; precs]
    vals[2:end,1] = 1:length(tests)
    vals[2:end,2] = correct_digits.([taylor_2f1(test..., order = order) for test ∈ tests], tru)
    for (i, prec) ∈ enumerate(precs)
        setprecision(prec)
        vals[2:end, i + 2] = correct_digits.([taylor_2f1((big.(test))..., order = order) for test ∈ tests], tru)
    end

    return vals
end

relative_error(val::Number, tru::Number)::Float64 = abs((val - tru) / tru)

function correct_digits(val::Number, tru::Number)
    relerr = relative_error(val, tru)

    if iszero(relerr)
        digs = 15
    elseif isinf(relerr) || isnan(relerr)
        digs = 0
    else
        digs = -round(Int,log10(relerr))
        if digs <= 0 
            digs = 0
        elseif digs > 16
            digs = 16
        end
    end

    return digs
end
