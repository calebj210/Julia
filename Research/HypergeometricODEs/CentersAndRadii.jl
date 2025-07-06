#=
#   Routines for computing centers and radii of transformations
#
# Author: Caleb Jacobs
# DLM: July 6, 2025
=#

function z_rc(a, b, c, maxr)
    R = min(abs(c / a / b), maxr)

    return (; r = R, c = 0, d = 1)
end

function zalt_rc(a, b, c, maxr)
    R = min(abs(c / (c - a) / (c - b)), maxr)

    return (; r = R, c = 0, d = 1)
end

function zoverzminus1a_rc(a, b, c, maxr)
    R = min(abs(c / a / (c - b)), maxr)

    return (; r = R / (1 - R^2), c = R^2 / (R^2 - 1), d = 1)
end

function zoverzminus1b_rc(a, b, c, maxr)
    R = min(abs(c / (c - a) / b), maxr)

    return (; r = R / (1 - R^2), c = R^2 / (R^2 - 1), d = 1)
end

function oneminusz_rc(a, b, c, maxr)
    R = min(
        abs((a + b - c + 1) / a / b),
        abs((c - a - b + 1) / (c - a) / (c - b)),
        maxr
    )

    return (; r = R, c = 1, d = 1)
end

function oneminusoneoverz_rc(a, b, c, maxr)
    R = min(
        abs((a + b - c + 1) / a / (a - c + 1)),
        abs((c - a - b + 1) / (c - a) / (1 - a)),
        maxr
    )

    return (; r = 1 / (1 - R^2), c = 1 / (1 - R^2), d = 1)
end

function oneoverz_rc(a, b, c, maxr)
    R = min(
        abs((a - b + 1) / a / (a - c + 1)),
        abs((b - a + 1) / b / (b - c + 1)),
        maxr
    )

    return (; r = 1 / R, c = 0, d = -1)
end

function oneoveroneminusz_rc(a, b, c, maxr)
    R = min(
        abs((a - b + 1) / a / (c - b)),
        abs((b - a + 1) / b / (c - a)),
        maxr
    )
    
    return (; r = 1 / R, c = 1, d = -1)
end

function get_rcs(a, b, c, z, maxr = 0.5)
    # Transformations for which z falls within the reliable radius of convergence
    rcs = (;
        z = z_rc(a, b, c, maxr),
        zalt = zalt_rc(a, b, c, maxr),
        zoverzminus1a = zoverzminus1a_rc(a, b, c, maxr),
        zoverzminus1b = zoverzminus1b_rc(a, b, c, maxr),
        oneminusz = oneminusz_rc(a, b, c, maxr),
        oneminusoneoverz = oneminusoneoverz_rc(a, b, c, maxr),
        oneoverz = oneoverz_rc(a, b, c, maxr),
        oneoveroneminusz = oneoveroneminusz_rc(a, b, c, maxr),
    )

    in_radii = Vector{Symbol}()
    for (k, rc) ∈ pairs(rcs)
        if rc.d == 1
            if abs(z - rc.c) < rc.r
                push!(in_radii, k)
            end
        else
            if abs(z - rc.c) > rc.r
                push!(in_radii, k)
            end
        end
    end

    return (rcs, in_radii)
end
