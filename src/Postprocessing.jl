module Postprocessing

using LinearAlgebra
using FFTW
using OffsetArrays
using ToeplitzMatrices

include("total_variation.jl")
include("fourier_pade.jl")

export total_variation, total_variation_denoising, total_variation_denoising!,
       group_sparse_total_variation_denoising, group_sparse_total_variation_denoising!
export fourier_pade

end # module
