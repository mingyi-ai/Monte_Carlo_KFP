using Random: Xoshiro, randn

# Every trajectory receives a separately seeded standard-library generator, so
# results do not depend on thread scheduling. Accelerator backends should use a
# device-appropriate counter-based generator and validate statistical agreement.
@inline function splitmix64(value::UInt64)
    z = value + 0x9e3779b97f4a7c15
    z = (z ⊻ (z >> 30)) * 0xbf58476d1ce4e5b9
    z = (z ⊻ (z >> 27)) * 0x94d049bb133111eb
    return z ⊻ (z >> 31)
end

@inline path_rng(seed::UInt64, trajectory::Int) =
    Xoshiro(splitmix64(seed + UInt64(trajectory)))

@inline normal!(rng::Xoshiro, ::Type{T}) where {T<:AbstractFloat} =
    randn(rng, T)
