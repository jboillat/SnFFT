
# Parameters:
#	N::Int
#	- the problem size
#	M::Int
#	- the module
#	PA::Array{Array{Int, 1}, 1}
#	- PA[i] is a Permutation of N
#	VA::Array{Bigint, 1}
#	- VA[i] is the Value associated with PA[i]
# Return Values:
#	SNF::Array{BigInt, 1}
#	- SNF[i] is the value associated with the permutation that permutation_index() maps to i
#	- this is the format for the SNF parameter of sn_fft()
# Notes:
#	- any permutation of N not represented in PA will be assigned a value of zero
function snf(N::Int, M::Int, PA::Array{Array{Int64, 1}, 1}, VA::Array{BigInt, 1})
        SNF = zeros(BigInt, factorial(N))
	for i = 1:length(PA)
		Permutation = PA[i]
		Index = permutation_index(Permutation)
		SNF[Index] = big(VA[i]) % M
	end
	return SNF
end
