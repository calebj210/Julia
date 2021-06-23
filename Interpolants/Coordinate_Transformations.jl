#=
  Methods to convert from one coordinate system to another

  Date last modified: 12-07-2020
  Author: Caleb Jacobs
=#

"""
    Convert Cartesiian data to polar data

    data is the data to converted
"""
function cartToPolar(data)
    polData = copy(data)

    polData[:,1] = atan.(data[:,2], data[:,1])

    for i = 1:size(data,1)
        polData[i,2] = norm(data[i,:])
    end

    return polData
end

"""
    Convert polar data to Cartesian data

    data is the data to be converted
"""
function polarToCart(data)
    cartData = hcat(cos.(data[:,1]), sin.(data[:,1]))

    cartData .*= data[:,2]

    return cartData
end

"""
    Convert Cartesian data to exponential data (x -> x, y -> sign(y)ℯˣ-1)

    data is the data to be converted
"""
function cartToExp(data)
    expData = hcat(data[:,1], log.(data[:,2] .+ 1))

    return expData
end

"""
    Convert exponential data to Cartesian data

    data is the data to be converted
"""
function expToCart(data)
    cartData = hcat(data[:,1], expm1.(data[:,2]))

    return cartData
end

"""
    Convert Cartesian data to sinusoidal data

    data is the data to be converted
"""
function cartToSin(data)
    cartData = hcat(data[:,1], data[:,2] - sin.(data[:,1]))

    return cartData
end

"""
    Converte sinusoidal data to Cartesian data

    data is the data to be converted
"""
function sinToCart(data)
    sinData = hcat(data[:,1], sin.(data[:,1]) + data[:,2])

    return sinData
end
