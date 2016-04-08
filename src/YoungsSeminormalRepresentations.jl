# Parameters
#	N::Int
#	- the problem size
#       M::Int
#       - the module
# Return Values
#	YSR::Array{Array{Array{SparseMatrixCSC, 1}, 1}, 1} (Young's Seminormal Representations)
#	- YSR[n][p][k] is Young's Seminormal Representation for the Adjacent Transposition (K, K + 1)
#         for the pth Partition of n
#	PT::Array{Array{Array{Int, 1}, 1}, 1} (Partition Tree)
#	- output1 from partition_tree()
function ysr(N::Int, M::Int)
	P, WI = partitions(N)
	PT = partition_tree(N, P, WI)
	YS = ys_symbols(N, P, PT)
	YSR = Array(Array{Array{SparseMatrixCSC, 1}, 1}, N)
	IYSR = Array(Array{Array{SparseMatrixCSC, 1}, 1}, N)
	for n = 1:N
		YSn = YS[n]
		Pn = P[n]
		Pn_L = length(Pn)
		YSRn = Array(Array{SparseMatrixCSC, 1}, Pn_L)
		for p = 1:Pn_L
			YSnp = YSn[p]
			Pnp = Pn[p]
			R = length(Pnp)
			YSRn[p] = YSR_p(n, M, R, YSnp)
		end
		YSR[n] = YSRn
	end
	return YSR, PT
end

# Parameter
#	N::Int
#	- the size of P
#       M::Int
#       - the module
#	R::Int
#	- the length of P
#	P::Array{Int, 1}
#	- P is a Partition of N
# Return Values:
#	YSRp::Array{SparseMatrixCSC, 1}
#	- YSRp[k] is the Young's Seminormal Representation of P for the Adjacent Transposistion (k, k + 1)
function YSR_p(N::Int, M::Int, R::Int, YSymbols::Array{Array{Int, 1}, 1})
	YSRp = Array(SparseMatrixCSC, N - 1)
	L = length(YSymbols)
	DA, IA = ys_information(N, R, YSymbols)
	for k = 1:(N - 1)
		YSRp[k] = YSR_pk(M, DA, IA, L, k)
	end
	return YSRp
end

#Parameters
#       M::Int
#       - the module
#	DA::Array{Int, 2}
#	- output1 from ys_information()
#	IA::Array{Int, 2}
#	- outpu2 from ys_information()
#	L::Int
#	- the number of Yamanouchi Symbols for the current Partition
#	K::Int
#	- represents the Adjacent Transposition (K, K + 1)
#Return Values
#	YSRpk::SparseMatrixCSC
#	- Young's Seminormal Representation of the Adjacent Transposition (K, K + 1) for the Partition p
function YSR_pk(M::Int, DA::Array{Int, 2}, IA::Array{Int, 2}, L::Int, k::Int)
	n = L
	for i = 1:L
		if IA[i,k] != 0
			n += 1
		end
	end
	colptr = Array(Int, L + 1)
	colptr[1] = 1
	rowval = Array(Int, n)
	nzval = Array(BigInt, n)
	cp = 1
	for i = 1:L
	        Dik = DA[i, k]
		D = invmod(big(Dik),M)
		I = IA[i,k]
		if I == 0
			colptr[i + 1] = colptr[i] + 1
			rowval[cp] = i
			nzval[cp] = big((M + sign(Dik)) % M) # [original SnFFT] nzval[cp] = D
			cp += 1
		else
			colptr[i + 1] = colptr[i] + 2
			if I > i
				rowval[cp] = i
				nzval[cp] = big(D)
				cp += 1
				rowval[cp] = I
				nzval[cp] = big(1) # [original SnFFT] nzval[cp] = sqrt(1 - D * D)
				cp += 1
			else
				rowval[cp] = I
        nzval[cp] = (M + (M + 1 - D * D % M) % M) % M # [original SnFFT] nzval[cp] = sqrt(1 - D * D)
				cp += 1
				rowval[cp] = i
				nzval[cp] = big(D)
				cp += 1
			end
		end
	end
	YSRpk = SparseMatrixCSC(L, L, colptr, rowval, nzval)
	return YSRpk
end
