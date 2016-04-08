
# Parameters:
#       M::Int
#	- the module
#	FFT1::Array{Array{BigInt, 2}, 1}
#	- the first Fourier transform
#	FFT2::Array{Array{BigInt, 2}, 1}
#	- the second Fourier transform
# Return Values:
#	Convolution::Array{Array{BigInt, 2}, 1}
#	- the convolution of FFT1 and FFT2
# Notes:
#	- FFT1 and FFT2 have to be Fourier transforms of functions over the same Sn
#	- However, the don't have to have the same degree of bandlimitedness
function sn_convolution(M::Int, FFT1::Array{Array{BigInt, 2}, 1}, FFT2::Array{Array{BigInt, 2}, 1})
	L = length(FFT1)
	Convolution = Array(Array{BigInt, 2}, L)
	while L != 0
		Convolution[L] = ((FFT1[L] % M) * (FFT2[L] % M)) % M 
		L -= 1
	end
	return Convolution
end

