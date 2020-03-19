using Postprocessing
using Test
using LinearAlgebra, Random

@testset "total variation denoising" begin
  for T in (Float32, Float64)
    x = range(zero(T), 2*one(T), length=1000)
    y_true = sinpi.(x)
    y_noisy = y_true .+ 0.1 .* randn(T, length(x))
    y_denoised = total_variation_denoising(y_noisy, 1)
    @test norm(y_denoised - y_true) < norm(y_noisy - y_true)
    @test_throws ArgumentError total_variation_denoising(y_noisy, -1)
    @test_throws DimensionMismatch total_variation_denoising!(similar(y_noisy, T, length(y_noisy)-1), y_noisy, 1)
  end
end

@testset "group sparse total variation denoising" begin
  for T in (Float32, Float64)
    x = range(zero(T), 2*one(T), length=1000)
    y_true = sinpi.(x)
    y_noisy = y_true .+ 0.1 .* randn(T, length(x))
    y_denoised = group_sparse_total_variation_denoising(y_noisy, 1, group_size=3)
    @test norm(y_denoised - y_true) < norm(y_noisy - y_true)
    @test_throws ArgumentError group_sparse_total_variation_denoising(y_noisy, -1, group_size=3)
    @test_throws DimensionMismatch group_sparse_total_variation_denoising!(similar(y_noisy, T, length(y_noisy)-1), y_noisy, 1)
  end
end
