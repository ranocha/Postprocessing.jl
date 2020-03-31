using Postprocessing
using Test
using LinearAlgebra, DelimitedFiles

@testset "Fourier Padé reconstruction" begin
  for T in (Float32, Float64)
    data = readdlm(joinpath(@__DIR__, "test_u_piecewise_constant.txt"), comments=true)
    x = T.(data[:, 1])
    u = T.(data[:, 2])
    u_true = @. ifelse(abs(x - 0.6) < 0.2, one(T), zero(T))

    u_reconstructed = fourier_pade(x, u, 50, 50)
    @test norm(u_reconstructed - u_true) < norm(u - u_true)
    @test total_variation(u_reconstructed) < total_variation(u)

    u_reconstructed = fourier_pade(x, u, 50, 45)
    @test norm(u_reconstructed - u_true) < norm(u - u_true)
    @test total_variation(u_reconstructed) < total_variation(u)

    # this choice of the degrees leads to increased oscillations
    # near the discontinuity but should still be computable
    u_reconstructed = fourier_pade(x, u, 50, 52)

    @test_throws ArgumentError fourier_pade(x, u[1:end-1], 50, 50)

    N = length(u)
    @test_throws ArgumentError fourier_pade(x, u[1:end-1], N-(N÷2)+1, N÷2)
  end
end
