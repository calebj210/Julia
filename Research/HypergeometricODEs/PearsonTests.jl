include("pFq.jl")

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

"""
    run_pearson_tests(prec; order = 1000)
Run 2F1 tests mentioned in Pearson et al using the Taylor method of `order` and rerun test with extended `prec`ision calculations.
"""
function run_pearson_tests(prec = 256; order = 1000, relative = false)
    tru = [mathematica_2f1(test...) for test ∈ tests]
    valdouble = [taylor_2f1(test..., order = order) for test ∈ tests]

    setprecision(prec)
    valextend = [taylor_2f1((big.(test))..., order = order) for test ∈ tests]

    vals = [valdouble valextend]

    if relative
        return [[1:length(tests)...] abs.((vals .- tru) ./ (abs.(tru) .+ eps()))]
    else
        return [[1:length(tests)...] abs.(vals .- tru)]
    end
end
