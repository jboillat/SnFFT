# Parameters:
#	Permutation::Array{Int, 1}
#	- a permutation
# Return Values:
#	ST::String
#	- the string representation of permutation
function permutation_string(Permutation::Array{Int, 1})
	ST = "[ "
	for i = 1:length(Permutation)
		ST = string(ST, Permutation[i], " ")
	end
	ST = string(ST, "] ")
	return ST
end

# Parameters:
#	Partition::Array{Int, 1}
#	- a partition
# Return Values:
#	ST::String
#	- the string representation of partition
function partition_string(Partition::Array{Int, 1})
	ST = "( "
	for i = 1:length(Partition)
		ST = string(ST, Partition[i], " ")
	end
	ST = string(ST, ") ")
	return ST
end

# Parameters:
#	YS::Array{Int, 1}
#	- a yamanouchi symbol
# Return Values:
#	None - prints YS
function print_ys(YS::Array{Int, 1})
	N = length(YS)
	numRows = maximum(YS)
	ST = Array(AbstractString, numRows)
	IA = ones(Int, numRows)
	for i = 1:numRows
		ST[i] = ""
	end
	for n = 1:N	
		row = YS[n]
		NI = IA[row]
		ST[row] = string(ST[row]," ", n)
		NI += 1
		IA[row] = NI
	end
	for i = 1:numRows
		println(ST[i])
	end
end
