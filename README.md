# Postprocessing

[![Build Status](https://github.com/ranocha/Postprocessing.jl/workflows/CI/badge.svg)](https://github.com/ranocha/Postprocessing.jl/actions)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/ranocha/Postprocessing.jl?svg=true)](https://ci.appveyor.com/project/ranocha/Postprocessing-jl)
[![Codecov](https://codecov.io/gh/ranocha/Postprocessing.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/ranocha/Postprocessing.jl)
[![Coveralls](https://coveralls.io/repos/github/ranocha/Postprocessing.jl/badge.svg?branch=master)](https://coveralls.io/github/ranocha/Postprocessing.jl?branch=master)
[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
[![GitHub commits since tagged version](https://img.shields.io/github/commits-since/ranocha/Postprocessing.jl/v0.1.0.svg?style=social&logo=github)](https://github.com/ranocha/Postprocessing.jl)
<!-- [![PkgEval](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/P/Postprocessing.svg)](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/report.html) -->

This is work in progress. Currently, the following functions are implemented.
- `total_variation_denoising(source::AbstractVector, λ::Number)`
  Compute the solution of the total variation regularized least square problem
  ```math
      \\min_x \\frac{1}{2} \\sum_{k} |y_k - x_k|^2 + \\lambda \\sum_{k} |x_{k+1} - x_{k}|
  ```
  using the explicit algorithm of
  [Condat (2013) A Direct Algorithm for 1-D Total Variation Denoising](https://doi.org/10.1109/LSP.2013.2278339).
  An inplace version `total_variation_denoising!(dest::AbstractVector, source::AbstractVector, λ::Number)`
  is also provided.
- `group_sparse_total_variation_denoising(y::AbstractVector, λ::Number; group_size::Integer=1, max_iter::Integer=100)`
  Compute `max_iter` iterations of the algorithm described by
  [Selesnick and Chen (2013) Total variation denoising with overlapping group sparsity](https://doi.org/10.1109/ICASSP.2013.6638755).
- `fourier_pade(u, degree_num, degree_den, num_output=length(u))`
  Compute the Fourier-Padé reconstruction of `u` with degrees
  `(degree_num, degree_den)` and evaluate it at `num_output`
  equispaced points,
  cf. [Driscoll and Fornberg (2001) A Padé-based algorithm for overcoming the Gibbs phenomenon](https://doi.org/10.1023/A:1016648530648).
