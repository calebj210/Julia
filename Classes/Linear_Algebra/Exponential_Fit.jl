using LinearAlgebra
using LaTeXStrings
using Plots

### Linear exponential data fitter for MA4330 homework 9
function fitExp(data)
    A = [ℯ.^(data[:,1]) ℯ.^(-data[:,1])]
    b = data[:,2]

    x = (A' * A) \ (A' * b)

    return x
end

function main()
    data = [0 3.1;
            1 1.4;
            2 1.0;
            3 2.2;
            4 5.2;
            5 15.0]

    x = fitExp(data)

    t = range(0, 5, length = 100)
    f = x[1]*(ℯ.^t) + x[2]*(ℯ.^(-t))

    plotA = plot(t,f,
                 label = L"f(x) = c_1 e^x + c_2 e^{-x}")
    scatter!(data[:,1],data[:,2],
             markerstrokalpha = 0,
             markerstrokewidth = 0,
             markersize = 3,
             legend = :topleft,
             label = "Data",
             title = "Least Squares Data Fit Curve",
             xlabel = "x",
             ylabel = "y",
             dpi = 300)

    png(plotA,"ExpFit.png")

    display(plotA)

    print("\nx1 = $(x[1])\nx2 = $(x[2])\n")
end

main()
