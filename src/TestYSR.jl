include("Counter.jl")
include("Element.jl")
include("Hook.jl")
include("Partitions.jl")
include("PartitionTree.jl")
include("PrintUtils.jl")
include("YamanouchiSymbols_BL.jl")
include("YamanouchiSymbols.jl")
include("YoungsSeminormalRepresentations.jl")
println(">  YSR test")
function testYSR()
	testYSR(6,997,[1,1,1,1,1,1],[2,3,4,5,6,1])
	testYSR(6,997,[2,1,1,1,1],[2,3,4,5,6,1])
	testYSR(6,997,[2,2,1,1],[2,3,4,5,6,1])
	testYSR(6,997,[2,2,2],[2,3,4,5,6,1])
	testYSR(6,997,[3,1,1,1],[2,3,4,5,6,1])
	testYSR(6,997,[3,2,1],[2,3,4,5,6,1])
	testYSR(6,997,[4,1,1],[2,3,4,5,6,1])
	testYSR(6,997,[4,2],[2,3,4,5,6,1])
	testYSR(6,997,[5,1],[2,3,4,5,6,1])
	testYSR(6,997,[6],[2,3,4,5,6,1])
	
end

# Parameters:
#	N::Int
#	- the problem size
#       M::Int64
#       - the module
#	partition::Array{Int, 1}
#	- a partition of N
#	permutation::Array{Int, 1}
#	- a permutation of N
function testYSR(N::Int, M::Int64, partition::Array{Int, 1}, permutation::Array{Int, 1})
	P, WI = partitions(N) #Creates the partitions of i for i is 1 through N
 	p = 1 #Find the index of the specified partitions  
	while P[N][p] != partition
	  p += 1
	end

	YS = ys_partition(N, partition) #Creates the Yamanouchi Symbols of the specified partition
	degree = length(YS)

	YSR, PT = ysr(N,M) #Creates Young's Seminormal Representations and the Pointer Tree for partition decomposistion
	YSRnp = YSR[N][p] #Get the Representation of each Adjacent Transposition of the specified partition

	RM = ysr_permutation(M, permutation, YSRnp) #Creates the Representation Matrix (and its inverse) of the specified permutation for the specified partition

	println("Finds the representation of a permutation for a particular partition")
	println("")
	println("Permutation: ")
	println(permutation_string(permutation))
	println("")
	println("Partition: ")
	println(partition_string(partition))
	println("")
	println("Degree: ")
	println(degree)
	println("")
	println("Standard Tableaux: ")
	for i = 1:degree #Converts the Yamanouchi Symbols into Standard Tableau and prints them
		print_ys(YS[i])
		println("")
	end
	println("Representation matrix of the permutation for this partition: ")
	print(full(RM))
	println()
	println()
end



testYSR()
