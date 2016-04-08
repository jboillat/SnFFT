# Parameters:
#	N::Int
#	- the problem size
#	M::Int
#	- the module
#	SNF1::Array{Int, 1}
#	- SNF1[i] is the value associated with the Permutation that permutation_index() maps to i
#	SNF2::Array{Int, 1}
#	- SNF2[i] is the value associated with the Permutation that permutation_index() maps to i
#	YSR::Array{Array{Array{SparseMatrixCSC, 1}, 1}, 1}
#	- output1 from ysr()
#	PT::Array{Array{Array{Int, 1}, 1}, 1}
#	- output2 from ysr()
# Return Values:
#	C::Array{BigInt, 1}
#	- C the convolution product SNF1 * SNF2
function convolution(N::Int, M::Int, SNF1::Array{BigInt, 1}, SNF2::Array{BigInt, 1}, YSR::Array{Array{Array{SparseMatrixCSC, 1}, 1}, 1}, PT::Array{Array{Array{Int, 1}, 1}, 1}) 
  FFT1 = sn_fft(N, M, SNF1, YSR, PT) 
  FFT2 = sn_fft(N, M, SNF2, YSR, PT)
  @assert length(FFT1) == length(FFT2)
  Convolution = sn_convolution(M, FFT1,FFT2)
  C = sn_ifft(N, M, Convolution, YSR, PT)
  return C
end

# Parameters:
#	N::Int
#	- the problem size
#	M::Int
#	- the module
#	SNF::Array{BigInt, 1}
#	- SNF[i] is the value associated with the Permutation that permutation_index() maps to i
#       E::Int
#       - the exponent
# 	YOR::Array{Array{Array{SparseMatrixCSC, 1}, 1}, 1}
#	- output1 from ysr()
#	PT::Array{Array{Array{Int, 1}, 1}, 1}
#	- output2 from ysr()
# Return Values:
#	SNF1::Array{BigInt, 1}
#	- ISNF the convolution product SNF1 * SNF2
function fastExp(N::Int, M::Int, SNF::Array{BigInt, 1}, E::Int, YSR::Array{Array{Array{SparseMatrixCSC, 1}, 1}, 1}, PT::Array{Array{Array{Int, 1}, 1}, 1})
  P = factorial(N)
  R = zeros(BigInt, P)
  t = zeros(BigInt, P)
  for i = 1:P 
    t[i] = SNF[i]
  end  
  Id = zeros(Int, N)
  for i = 1:N
    Id[i] = i
  end
  R[permutation_index(Id)] = big(1)
  while  E > 0
    if E % 2 == 1
      u = convolution(N, M, R, t, YSR, PT)
      E -= 1
      for i = 1:P
        R[i] = u[i]
      end
    end
    u = convolution(N, M, t, t, YSR, PT)
    for i = 1:P
      t[i] = u[i]
    end
    E /= 2
  end
  return R
end
