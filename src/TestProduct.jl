include("Counter.jl")
include("Element.jl")
include("Hook.jl")
include("Partitions.jl")
include("PartitionTree.jl")
include("PrintUtils.jl")
include("YamanouchiSymbols_BL.jl")
include("YamanouchiSymbols.jl")
include("YoungsSeminormalRepresentations.jl")
println("> Product test")

#If a and b are permutations and c = a * b, demonstrates that the representation of c is the representation of a multiplied by the representation of b
function testProduct()
	testProduct(5,32416190071,[3,2],[1,2,4,3,5],[2,3,1,4,5])
end

# Parameters:
#	N::Int
#	- the problem size
#	M::Int
#	- the module
#	partition::Array{Int, 1}
#	- a partition of N
#	p1::Array{Int, 1}
#	- the first permutation of N
#	p2::Array{Int, 1}
#	- the second permutation of N
function testProduct(N::Int, M::Int, partition::Array{Int, 1}, p1::Array{Int, 1}, p2::Array{Int, 1})

	P, WI = partitions(N) #Creates the partitions of i for i is 1 through N
	p = 1 #Find the index of the specified partitions
	while P[N][p] != partition
		p += 1
	end

	pm = sn_multiply(p1, p2) #Find the product of the two permutations

	YSR, PT = ysr(N,M) #Creates Young's Orthogonal Representations and the Pointer Tree for partition decomposistion
	YSRnp = YSR[N][p] #Get the Representation of each Adjacent Transposition of the specified partition

	RM1, IRM1 = ysr_permutation(M, p1, YSRnp) #Find the representation of the first permutation
	RM2, IRM2 = ysr_permutation(M, p2, YSRnp) #Find the representation of the second permutation
	RMm, IRMm = ysr_permutation(M, pm, YSRnp) #Find the representation of the product of the permutations
	error = RM1 * RM2 - RMm #Find the error in the process

	println("If a and b are permutations and c = a * b, demonstrates that the representation of c is the representation of a multiplied by the representation of b")
	println("")
	ST1 = permutation_string(p1)
	println("Representation of ", ST1)
	println(RM1)
	println("")
	ST2 = permutation_string(p2)
	println("Representation of ", ST2)
	println(RM2)
	println("")
	STm = permutation_string(pm)
	println(ST1," * ", ST2, " is ", STm)
	println("")
	println("Representation of ", STm)
	println(RMm)
	println("")
	println("Product of ", ST1, " with ", ST2)
	println((RM1 * RM2) % M)
end

testProduct()
