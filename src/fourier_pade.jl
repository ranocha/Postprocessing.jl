"""
    fourier_pade(u, degree_num, degree_den, num_output=length(u))

Compute the Fourier-Padé reconstruction of `u` with degrees
`(degree_num, degree_den)` and evaluate it at `num_output`
equispaced points, cf.
Driscoll and Fornberg (2001) A Padé-based algorithm for overcoming the Gibbs phenomenon,
doi: 10.1023/A:1016648530648.
"""
function fourier_pade(u, degree_num, degree_den, num_output=length(u))
  N = degree_num + degree_den
  modes = length(u) ÷ 2 + 1
  if modes <= N
    throw(ArgumentError("Cannot perform a Fourier-Padé reconstruction with degrees ($degree_num, $degree_den) given $(length(u)) real coefficients corresponding to $modes complex modes."))
  end

  uhat = rfft(u)
  uhat ./= length(u)

  # denominator
  col_den = uhat[degree_num+2:N+1]
  row_den = zeros(eltype(col_den), degree_den+1)
  row_den[1:min(degree_num+2,degree_den+1)] = uhat[degree_num+2:-1:max(1, degree_num+2-degree_den)]
  Z = nullspace(Matrix(Toeplitz(col_den, row_den)))
  den_p = zeros(eltype(uhat), num_output)
  den_p[1:degree_den+1] = Z[:, end]
  den_p ./= den_p[findfirst(!iszero, den_p)]

  # numerator
  col_num = uhat[1:degree_num+1]
  col_num[1] /= 2
  row_num = zeros(eltype(col_num), degree_den+1)
  row_num[1] = col_num[1]
  A = Toeplitz(col_num, row_num)
  num_p = zeros(eltype(uhat), num_output)
  mul!(view(num_p, 1:degree_num+1), A, view(den_p, 1:degree_den+1))

  # evaluate Fourier-Padé approximation
  bfft_plan = plan_bfft(num_p)
  2 .* real.((bfft_plan * num_p) ./ (bfft_plan * den_p))
end


function reference(::typeof(fourier_pade))
"""
@article{driscoll2001pade,
  title={A {P}ad{\'e}-based algorithm for overcoming the {G}ibbs phenomenon},
  author={Driscoll, Tobin A and Fornberg, Bengt},
  journal={Numerical Algorithms},
  volume={26},
  number={1},
  pages={77--92},
  year={2001},
  publisher={Springer},
  doi={10.1023/A:1016648530648}
}
"""
end
