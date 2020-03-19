# Postprocessing

[![Build Status](https://travis-ci.com/ranocha/Postprocessing.jl.svg?branch=master)](https://travis-ci.com/ranocha/Postprocessing.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/ranocha/Postprocessing.jl?svg=true)](https://ci.appveyor.com/project/ranocha/Postprocessing-jl)
[![Codecov](https://codecov.io/gh/ranocha/Postprocessing.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/ranocha/Postprocessing.jl)
[![Coveralls](https://coveralls.io/repos/github/ranocha/Postprocessing.jl/badge.svg?branch=master)](https://coveralls.io/github/ranocha/Postprocessing.jl?branch=master)

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
