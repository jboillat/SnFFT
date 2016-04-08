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

println(">  Diffie-Hellman test")


# Demonstrates DH on Sn
# Parameters:
#	N::Int
#	- the problem size
#	M::Int
#	- the module
# Return
#	Nothing
function testDH(N::Int, M::Int)  
	YSR, PT = ysr(N, M) #Creates Young's Seminormal Representations and the Pointer Tree for Partition decomposistion
	alpha = zeros(BigInt, factorial(N))
	sum = 0
        for i = 2:factorial(N)
          alpha[i] = big(rand(0:M-1))
          sum += alpha[i]
        end
        alpha[1] = -sum
        alpha = (alpha % M + M) % M
	println("")
	println()
	println("alpha ", alpha)
	println()
#	s = 1327329843984343234
#	t = 1988328326632633323
	s = 132
	t = 19
	t1 = time()
	publicAlice = fastExp(N, M, alpha, s, YSR, PT)
	publicBob = fastExp(N, M, alpha, t, YSR, PT)
        privateAlice = fastExp(N, M, publicBob, s, YSR, PT)
	privateBob = fastExp(N, M, publicAlice, t, YSR, PT)
        t2 = time()
        println("Alice ", privateAlice)
        println("Bob ", privateBob)
	@assert maximum(abs(privateAlice -privateBob)) == 0
	println()
	println("elapsed time: ",t2-t1)
	sum = 0
end

p = 32416190071
#testDH(7,32416190071) # 154.66 s
#testDH(6,32416190071) # 8.55 s
#testDH(5,32416190071) # 0.83 s
#testDH(4,961748941) # 0.83 s
#testDH(3,32452867) # 0.11 s
#testDH(2,32452867) # 0.05 s

