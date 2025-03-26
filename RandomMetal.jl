# RandomMetal.jl
module RandomMetal

export xorshift32, xorshift32_float, box_muller

"""
    xorshift32(seed::UInt32) -> UInt32

A simple and efficient GPU-friendly XORSHIFT32 random number generator that
produces a pseudo-random 32-bit unsigned integer in the range [0, 2^32-1] from a given seed.
"""
function xorshift32(seed::UInt32)::UInt32
    seed ⊻= (seed << 13)
    seed ⊻= (seed >> 17)
    seed ⊻= (seed << 5)
    return seed
end

"""
    xorshift32_float(seed::UInt32) -> Float32

Converts the output of `xorshift32` to a Float32 number in (0,1).
It ensures a nonzero value by taking the maximum of the result and `1.0f-16`.
"""
function xorshift32_float(seed::UInt32)::Float32
    value = Float32(xorshift32(seed)) * 2.3283064f-10  # Scale to [0,1)
    return max(value, 1.0f-16)  # Ensure it's in (0,1)
end

"""
    box_muller(u1::Float32, u2::Float32) -> Float32

Generates a normally distributed random number from two uniformly distributed
Float32 numbers using the Box-Muller transform.
"""
function box_muller(u1::Float32, u2::Float32)
    r = sqrt(-2.0f0 * log(u1))
    theta = 2.0f0 * Float32(pi) * u2
    return r * cos(theta)
end

end  # module RandomMetal
