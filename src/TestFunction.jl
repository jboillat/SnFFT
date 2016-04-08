using Integers
include("Function.jl")
include("Element.jl")
include("PrintUtils.jl")

println("> Test Function")
#Demonstrates how to put the values of a function on Sn in the order used to compute the FFT
function testFunction(M::Int)

	PA = Array(Array{Int, 1}, 24) #Create an array to hold all of the permutations of 1,2,3,4
	PA[1] = [1,2,3,4]
	PA[2] = [1,2,4,3]
	PA[3] = [1,3,2,4]
	PA[4] = [1,3,4,2]
	PA[5] = [1,4,2,3]
	PA[6] = [1,4,3,2]
	PA[7] = [2,1,3,4]
	PA[8] = [2,1,4,3]
	PA[9] = [2,3,1,4]
	PA[10] = [2,3,4,1]
	PA[11] = [2,4,1,3]
	PA[12] = [2,4,3,1]
	PA[13] = [3,1,2,4]
	PA[14] = [3,1,4,2]
	PA[15] = [3,2,1,4]
	PA[16] = [3,2,4,1]
	PA[17] = [3,4,1,2]
	PA[18] = [3,4,2,1]
	PA[19] = [4,1,2,3]
	PA[20] = [4,1,3,2]
	PA[21] = [4,2,1,3]
	PA[22] = [4,2,3,1]
	PA[23] = [4,3,1,2]
	PA[24] = [4,3,2,1]

        VA = Array(BigInt, 24) #Create the values of the function on Sn
        for i = 1:24
          VA[i] = big(i)
        end
        VA[1] = parse(BigInt,"487989853985398530985385305303053050535385387538") % M
	SNF = snf(4,M,PA,VA) #Put the values into the order used to compute the FFT

	println("Demonstrates how to put the values of a function on Sn in the order used to compute the FFT")
	println("")
	println("Start with any ordering: ")
	println("Permutation   Function value for that Permutation")
	for i = 1:24
		ST = permutation_string(PA[i])
		ST = string(ST, "  ", VA[i])
		println(ST)
	end
	println("")
	println("Is sorted into this order for SnFFT to use: ")
	println("Permutation   Function value for that Permutation")
	for i = 1:24
		ni = SNF[i]
		ST = permutation_string(PA[i])
		ST = string(ST, "  ", SNF[i])
		println(ST)
	end
end

testFunction(997)