#=
#   Routines for initializing the Taylor method for 2F1
#
# Author: Caleb Jacobs
# DLM: July 23, 2025
=#

include("Transformations.jl")
include("CentersAndRadii.jl")

function intersectionpoint(z, rc)
    r = rc.r; c = rc.c

    if c == 0
        return r * sign(z)
    else
        proj = rc.c * real(z) / conj(z)
        if rc.d == -1 || abs2(proj - c) < r^2
            shift = sqrt(abs2(z) * r^2 - (imag(z) * c)^2) / conj(z)
            if  abs2(proj + shift - z) < abs2(proj - shift - z)
                return proj + shift
            else
                return proj - shift
            end
        else
            return Inf + Inf * im
        end
    end
end

function init(a, b, c, z, trans, rc)
    zs = NamedTuple{trans}([intersectionpoint(z,rc[tran]) for tran ∈ trans])

    tran = last(findmin(x -> abs(z - x), zs))
    z0 = zs[tran]
    fn = transformations[tran](a, b, c, z0)

    return (z0, fn)
end

function inward_init(a, b, c, z, rc)
    if 0 < real(z) && real(z) < 1
        trans = (:oneoverz, :oneminusoneoverz, :oneminusz)
        # trans = (:oneoverz, :oneoveroneminusz, :oneminusoneoverz, :oneminusz)
    else
        trans = (:oneoverz, :oneoveroneminusz)
    end

    return init(a, b, c, z, trans, rc)
end

function outward_init(a, b, c, z, rc)
    if real(z) > 1
        trans = (
            :z,
            :zalt,

            :zoverzminus1a,
            :zoverzminus1b,

            :oneminusz,
            :oneminusoneoverz,
        )
    else
        trans = (
            :z,
            :zalt,

            :zoverzminus1a,
            :zoverzminus1b,
        )
    end

    return init(a, b, c, z, trans, rc)
end

function initialize(a, b, c, z, maxr = 0.5)
    rcs, in_radii = get_rcs(a, b, c, z, maxr)

    if !isempty(in_radii)
        fn = transformations[last(findmax(x -> x.r, rcs[in_radii]))](a, b, c, z)

        return (z, fn)
    end

    if sign(real(a)) == sign(real(b)) && (real(a) < 0 || real(c) < 0)
        (z0, fn) = inward_init(a, b, c, z, rcs)
    else
        (z0, fn) = outward_init(a, b, c, z, rcs)
    end

    return (z0, fn)
end
