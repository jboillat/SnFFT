include("Counter.jl")
include("Convolution.jl")
include("Crypto.jl")
include("Element.jl")
include("FFT.jl")
include("Hook.jl")
include("IFFT.jl")
include("Partitions.jl")
include("PartitionTree.jl")
include("PrintUtils.jl")
include("YamanouchiSymbols_BL.jl")
include("YamanouchiSymbols.jl")
include("YoungsSeminormalRepresentations.jl")

println(">  Convolution test")


# Demonstrates fast convolution and standard convolution deliver the same result
# Parameters:
#	N::Int
#	- the problem size
#	M::Int
#	- the module
# Return
#	Nothing
function testConvolution(N::Int, M::Int)  
	YSR, PT = ysr(N, M) #Creates Young's Seminormal Representations and the Pointer Tree for Partition decomposistion
	alpha = zeros(BigInt, factorial(N))
	beta = zeros(BigInt, factorial(N))
	c2 = zeros(BigInt, factorial(N))
        for i = 1:factorial(N)
          alpha[i] = big((rand(Int) % 10 + 10) % 10)
          beta[i] = big((rand(Int) % 10 + 10) % 10)
        end
	# println("")
	# println()
	# println("alpha ", alpha)
	# println()
	# println("beta ", beta)
	# println()
	t1 = time()
	c1 = convolution(N, M, alpha, beta, YSR, PT)
	t2 = time()
	println("fast convolution ", t2-t1)
	println()  
	# println("c1 ", c1)
	# println()
	# standard slow convolution
        t1 = time()
	for i = 1:factorial(N)
	  pi = index_permutation(N,i)
	  for j = 1:factorial(N)
	    pj = index_permutation(N,j)
	    c2[i] = (c2[i] + (alpha[permutation_index(pj)] * beta[permutation_index(sn_multiply(sn_inverse(pj),pi))]) % M) % M
	  end
        end
        t2 = time()
	# println("c2 ", c2)
	# println()
	println("standard convolution ", t2-t1)
	println()  
	@assert maximum(abs(c1-c2)) == 0
end

testConvolution(7,997)
#testConvolution(9,997)
