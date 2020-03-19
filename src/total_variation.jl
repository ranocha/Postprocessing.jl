"""
    total_variation_denoising(y::AbstractVector, λ::Number)

Compute the solution of the total variation regularized least square problem
```math
    \\min_x \\frac{1}{2} \\sum_{k} |y_k - x_k|^2 + \\lambda \\sum_{k} |x_{k+1} - x_{k}|
```
using the explicit algorithm of
Condat (2013) A Direct Algorithm for 1-D Total Variation Denoising,
doi: 10.1109/LSP.2013.2278339.
"""
function total_variation_denoising(source::AbstractVector, λ::Number)
  T = promote_type(eltype(source), typeof(λ))
  dest = similar(source, T)
  total_variation_denoising!(dest, source, convert(T, λ))
end

"""
    total_variation_denoising!(dest::AbstractVector, source::AbstractVector, λ::Number)

In-place version of `total_variation_denoising(source, λ)` storing the result in `dest`.
"""
function total_variation_denoising!(dest::AbstractVector, source::AbstractVector, λ::Number)
  if λ < 0
    throw(ArgumentError("The regularization parameter λ must be non-negative!"))
  end

  idx = eachindex(dest, source)
  k = k₀ = k₋ = k₊ = first(idx)
  vmin = source[k] - λ
  vmax = source[k] + λ
  umin = λ
  umax = -λ

  if length(idx) == 1
    dest .= source
    return dest
  end

  @inbounds while true
    if k == length(idx)
      dest[idx[k]] = vmin + umin
      break
    end
    nextval = source[idx[k+1]]
    if nextval + umin < vmin - λ
      @simd for ki in k₀:k₋
        dest[idx[ki]] = vmin
      end
      k = k₀ = k₋ = k₊ = k₋ + 1
      thisval = source[idx[k]]
      vmin = thisval
      vmax = thisval + 2λ
      umin = λ
      umax = -λ
    elseif nextval + umax > vmax + λ
      @simd for ki in k₀:k₊
        dest[idx[ki]] = vmax
      end
      k = k₀ = k₋ = k₊ = k₊ + 1
      thisval = source[idx[k]]
      vmin = thisval - 2λ
      vmax = thisval
      umin = λ
      umax = -λ
    else
      k += 1
      umin += nextval - vmin
      umax += nextval - vmax
      if umin >= λ
        vmin += (umin - λ) / (k - k₀ + 1)
        umin = λ
        k₋ = k
      end
      if umax <= -λ
        vmax += (umax + λ) / (k - k₀ + 1)
        umax = -λ
        k₊ = k
      end
    end

    if k < length(idx)
      continue
    elseif umin < 0
      @simd for ki in k₀:k₋
        dest[idx[ki]] = vmin
      end
      k = k₀ = k₋ = k₋ + 1
      thisval = source[idx[k]]
      vmin = thisval
      umin = λ
      umax = thisval + λ - vmax
    elseif umax > 0
      @simd for ki in k₀:k₊
        dest[idx[ki]] = vmax
      end
      k = k₀ = k₊ = k₊ + 1
      thisval = source[idx[k]]
      vmax = thisval
      umin = thisval - λ - vmin
      umax = -λ
    else
      @simd for ki in k₀:length(idx)
        dest[idx[ki]] = vmin + umin / (k - k₀ + 1)
      end
      break
    end
  end

  dest
end
