Base.:(:)(start::Complex, step::Real, stop::Complex) = 
    StepRangeLen(
        step > 0 ? start : stop, 
        step * sign(stop - start), 
        floor(Int64, abs((stop - start) / step) + 1)
    )
Base.:(:)(start::A, step::Real, stop::C) where {A<:Complex, C<:Real} = Base.:(:)(convert(promote_type(A, C), start), step, convert(promote_type(A, C), stop))
Base.:(:)(start::A, step::Real, stop::C) where {A<:Real, C<:Complex} = Base.:(:)(convert(promote_type(A, C), start), step, convert(promote_type(A, C), stop))

Base.range(start::Complex, stop::Complex, length::Int64) = StepRangeLen(start, (stop - start) / (length - 1), length)
Base.range(start::A, stop::B, length::Int64) where {A<:Complex, B<:Real} = StepRangeLen(start, (stop - start) / (length - 1), length)
Base.range(start::A, stop::B, length::Int64) where {A<:Real, B<:Complex} = StepRangeLen(start, (stop - start) / (length - 1), length)
