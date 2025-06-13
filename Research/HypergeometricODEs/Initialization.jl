#=
#   Routines for initializing the Taylor method for 2F1
#
# Author: Caleb Jacobs
# DLM: June 5, 2025
=#

include("Transformations.jl")

function get_radii(a, b, c, z, maxr = 0.5)
    # Reliable transformation radii
    r = (;
        z = abs(c / a / b),
        zalt = abs(c / (c - a) / (c - b)),

        zoverzminus1a = abs(c / a / (c - b)),
        zoverzminus1b = abs(c / (c - a) / b),

        oneminusz = min(
            abs((a + b - c + 1) / a / b),
            abs((c - a - b + 1) / (c - a) / (c - b)),
        ),
        oneminusoneoverz = min(
            abs((a + b - c + 1) / a / (a - c + 1)),
            abs((c - a - b + 1) / (c - a) / (1 - a)),
        ),

        oneoverz = min(
            abs((a - b + 1) / a / (a - c + 1)),
            abs((b - a + 1) / b / (b - c + 1)),
        ),
        oneoveroneminusz = min(
            abs((a - b + 1) / a / (c - b)),
            abs((b - a + 1) / b / (c - a)),
        ),
    )

    # Transformations for which z falls within the reliable radius of convergence
    in_radii = findall((;
        z = abs(z) < min(maxr, r.z),
        zalt = abs(z) < min(maxr, r.zalt),

        zoverzminus1a = abs(z / (z - 1)) < min(maxr, r.zoverzminus1a),
        zoverzminus1b = abs(z / (z - 1)) < min(maxr, r.zoverzminus1b),

        oneminusz = abs(1 - z) < min(maxr, r.oneminusz),
        oneminusoneoverz = abs(1 - 1 / z) < min(maxr, r.oneminusoneoverz),

        oneoverz = abs(1 / z) < min(maxr, r.oneoverz),
        oneoveroneminusz = abs(1 / (1 - z)) < min(maxr, r.oneoveroneminusz),
    ))

    return (;r, in_radii)
end

function inward_init(a, b, c, z, r, maxr)
    if r.oneoverz >= r.oneoveroneminusz
        # 1 / z
        z0 = sign(z) / min(r.oneoverz, maxr)
        fn = oneoverz_2f1(a, b, c, z0)
    else
        # 1 / (1 - z)
        z0 = 1 + sign(z - 1) / min(r.oneoveroneminusz, maxr)
        fn = oneoveroneminusz_2f1(a, b, c, z0)
    end

    return (z0, fn)
end

function outward_init(a, b, c, z, r, maxr)
    if real(z) < 1
        if r.z >= r.zalt
            # z
            z0 = min(r.z, maxr) * sign(z)
            fn = maclaurin_2f1(a, b, c, z0)
        else
            # zalt
            z0 = min(r.zalt, maxr) * sign(z)
            fn = zalt_2f1(a, b, c, z0)
        end
    else
        if r.oneminusz >= r.oneminusoneoverz
            # 1 - z
            z0 = 1 + min(r.oneminusz, maxr) * sign(z - 1)
            fn = oneminusz_2f1(a, b, c, z0)
        else
            # 1 - 1/z
            R = min(r.oneminusoneoverz, maxr)               # Radius of disk of convergence
            center = 1 / (1 - R^2)                          # Center of disk of convergence
            z0 = center + R * sign(z - center)
            fn = oneminusoneoverz_2f1(a, b, c, z0)
        end
    end

    return (z0, fn)
end

function initialize(a, b, c, z, maxr = 0.5)
    radii, in_radii = get_radii(a, b, c, z, maxr)

    if !isempty(in_radii)
        fn = transformations[last(findmax(radii[in_radii]))](a, b, c, z)

        return (z, fn)
    end

    if sign(real(a)) == sign(real(b))
        if real(a) < 0
            (z0, fn) = inward_init(a, b, c, z, radii, maxr)
        else
            if real(c) > 0
                (z0, fn) = outward_init(a, b, c, z, radii, maxr)
            else
                (z0, fn) = inward_init(a, b, c, z, radii, maxr)
            end
        end
    else
        (z0, fn) = outward_init(a, b, c, z, radii, maxr)
    end

    return (z0, fn)
end
