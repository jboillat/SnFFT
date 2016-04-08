include("Counter.jl")
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

println(">  FFT test")
#Demonstrates how to compute the FFT, the IFFT, and shows that the IFFT recovers the initial function
function testFFT()
	testFFT(3,32416190071)
	#testFFT(4,997)
end

# Parameters:
#	N::Int
#	- the problem size
#       M::Int
#       - the module
function testFFT(N::Int, M::Int)

	YSR, PT = ysr(N,M) #Creates Young's Seminormal Representations and the Pointer Tree for Partition decomposistion
	SNF = zeros(BigInt, factorial(N))
        for i = 1:factorial(N)
          SNF[i] = Int(i)
        end
	FFT = sn_fft(N, M, SNF, YSR, PT) #Find the FFT
	iSNF = sn_ifft(N, M, FFT, YSR, PT) #Find the inverse FFT
	dif = SNF - iSNF #Find the error in the process

	P, WI = partitions(N) #Create the Partitions of 1:N
	Pn = P[N] #Get the partitions that act as the labels for the irreducibles


	println("Demonstrates how to compute the FFT, the IFFT, and shows that the IFFT recovers the initial function")
	println("")
	println("A random function on Sn: ")
	for i = 1:length(SNF)
		println(SNF[i])
	end
	println("")
	println("The Fourier Transform: ")
	for i = 1:length(FFT)
		println(partition_string(Pn[i]))
		println(FFT[i])
		println("")
	end
	println("The Inverse Fourier Transform: ")
	for i = 1:length(iSNF)
		println(SNF[i], "  ", iSNF[i] % M)
	end
	println("Maximum error in recovering the original function on Sn: ")
	println(maximum(abs(dif)))
end

testFFT()