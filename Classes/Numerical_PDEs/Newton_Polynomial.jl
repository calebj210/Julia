### 3rd Order Newton Polynomial in Julia ###

# General Newton polynomial weights calculator using divide differences
function compWeights(data)
    # Number of points
    n = size(data,2);

    # Copying data
    x = data[1,:];
    w = deepcopy(data[2,:]);

    # Computing divided differences
    for i ∈ 2:n
        for j ∈ 1:i-1
            w[i] = (w[j] -w[i])/(x[j] - x[i]);
        end
    end
    
    return w
end

# Main function for computing
function compPolynomial()
    # Given Data
    data = [-0.5 0 0.5 1;
            1 1.2 -0.2 0];
    
    # Compute weights
    w = compWeights(data);

    # Format results
    str = string("Output:\n\n",
                 "[x1,x2,x3,x4] = ", data[1,:], '\n',
                 "[a,b,c,d] = ", w, "\n\n",
                 "Final answer:\n",
                 "p(x) = a + b(x-x1) + c(x-x1)(x-x2) + d(x-x1)(x-x2)(x-x3)")

    # Store results
    io = open("tmp2", "w+");
    write(io,str);
    close(io);

    return w
end

compPolynomial()
