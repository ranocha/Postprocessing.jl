"""
    fourier_pade(x, u, degree_num, degree_den)

Compute the Fourier-Padé reconstruction of `u` with degrees
`(degree_num, degree_den)` and evaluate it at the points `x`, cf.
Driscoll and Fornberg (2001) A Padé-based algorithm for overcoming the Gibbs phenomenon,
doi: 10.1023/A:1016648530648.
"""
function fourier_pade(x, u, degree_num, degree_den)
  N0 = length(u)
  if iseven(N0)
    throw(ArgumentError("u has $N0 coeffficients but is assumed to have an odd number of coefficients."))
  end
  N = degree_num + degree_den
  if N0 <= N
    throw(ArgumentError("Cannot perform a Fourier-Padé reconstruction with degrees $degree_num and $degree_den given $N0 coefficients."))
  end

  uh = fft(u) / N0
  uhat = uh[1:N+1]

  # denominator
  col_den = uhat[degree_num+2:N+1]
  row_den = uhat[degree_num+2:-1:max(1, degree_num+2-degree_den)]
  if degree_den > degree_num
    row_den = vcat(row_den, fill(0, degree_den-degree_num-1))
  end
  Z = nullspace(Matrix(Toeplitz(col_den, row_den)))
  den_p = Z[:, end]
  den_p ./= den_p[findfirst(!iszero, den_p)]
  den_m = conj.(den_p)

  # numerator
  col_num = uhat[1:degree_num+1]
  col_num[1] /= 2
  row_num = zeros(eltype(col_num), degree_den+1)
  row_num[1] = col_num[1]
  A = Toeplitz(col_num, row_num)
  num_p = A * den_p
  num_m = conj.(num_p)

  # evaluate Fourier-Padé approximation
  T = real(eltype(uhat))
  n = length(x)
  xx = range(zero(T), 2*one(T)*π*(n-1)/n, length=n)
  x_p = @. exp(im * xx)
  x_m = @. exp(-im * xx)

  real.(evalpoly.(x_p, Ref(num_p)) ./ evalpoly.(x_p, Ref(den_p)) .+ evalpoly.(x_m, Ref(num_m)) ./ evalpoly.(x_m, Ref(den_m)))
end
