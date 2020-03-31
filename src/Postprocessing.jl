module Postprocessing

using LinearAlgebra
using FFTW
using OffsetArrays
using ToeplitzMatrices

"""
    reference(algorithm)

Returns a citable reference for `algorithm`, e.g. a bibtex entry
for a scientific article or a book.
"""
function reference end

include("total_variation.jl")
include("fourier_pade.jl")

export reference
export total_variation, total_variation_denoising, total_variation_denoising!,
       group_sparse_total_variation_denoising, group_sparse_total_variation_denoising!
export fourier_pade

end # module
