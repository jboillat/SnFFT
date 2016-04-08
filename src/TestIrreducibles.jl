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

println(">  Irreducible test")
#Demonstrates how to compute the irreducible representations
function testIrreducibles()
	testIrreducibles(4,997)
end

# Parameters:
#	N::Int
#	- the problem size
#       M::Int
#       - the module
function testIrreducibles(N::Int, M::Int)

	YSR, PT = ysr(N,M) #Creates Young's Seminormal Representations and the Pointer Tree for Partition decomposistion
	P, WI = partitions(N) #Create the Partitions of 1:N
	Pn = P[N] #Get the partitions that act as the labels for the irreducibles
	println("Demonstrates how to compute the irreducibles representations")
	println("")
	println("The irreducible representations: ")
	for i = 1:length(YSR)
		println(partition_string(Pn[i]))
		println()
		for j = 1:length(YSR[i]) 
		  for k = 1:length(YSR[i][j])
		    println(full(YSR[i][j][k]))
		    println()
                  end
                end
		println("")
	end
end

testIrreducibles()