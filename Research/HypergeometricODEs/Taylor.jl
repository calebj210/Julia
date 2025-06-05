#=
#   Taylor method routines for 2F1
#
# Author: Caleb Jacobs
# DLM: June 5, 2025
=#

function get_direction(z0, z, f, df, branch = false)
    # if branch && real(z) > 1 && real(z0) < 1 
    #     h_str = angle(1 + 2sgn(imag(z)) * im - z0)  # Go to 1 ± i first when navigating the branch point
    # else
    #     h_str = angle(z - z0)                       # Straight path direction
    #     branch = false
    # end

    # h_str = angle(z - z0)                       # Straight path direction
    # h_arg = angle(im * f / df)                  # Constant phase direction
    # if h_arg < 0
    #     h_arg += π
    # end
    #
    # h_dif = argmin(h -> abs2(h), h_arg - h_str .+ (-2:1) * π)
    #
    # return (cis(h_str + .25h_dif), branch)
    return (sign(z - z0), false)
end

function recursive_2f1(a, b, c, z0, f0, h, N; tol = eps())
    (c0, c1) = f0

    # 10 flop optimization for 3-term recurrence
    A0 = -a * b; A1 = (1 - a - b); A2 = -2
    B0 = B1 = (c - (1 + a + b) * z0); B2 = (2 - 4z0)
    C0 = C1 = C2 = 2z0 * (z0 - 1)
    
    hn = h
    S = c0 + c1 * h
    dS = c1
    crit_flag = false
    for n = 2:N
        # Compute next coefficient
        coeff = (A0 * c0 + B0 * c1) / C0

        # criteria = abs(val)
        # if criteria <= eps(abs(S))
        dS += coeff * hn * n
        hn *= h
        criteria = abs(coeff * hn / S)
        if criteria <= tol / 2
            if crit_flag
                break
            else
                crit_flag = true
            end
        elseif isnan(criteria) || isinf(criteria)
            # @warn "Method diverged"
            break
        elseif crit_flag
            crit_flag = false
        end
        S += coeff * hn

        c0 = c1
        c1 = coeff
        # push!(p.coeffs, coeff)

        # Update recurrence values
        A1 += A2; A0 += A1
        B1 += B2; B0 += B1
        C1 += C2; C0 += C1
    end 

    # return p
    return [S, dS]
end
