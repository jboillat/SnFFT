# Parameters:
#	N::Int
#	- the problem size
#	M::Int
#	- the module
#	SNF::Array{BigInt, 1}
#	- SNF[i] is the value associated with the Permutation that permutation_index() maps to i
#	YSR::Array{Array{Array{SparseMatrixCSC, 1}, 1}, 1}
#	- output1 from ysr()
#	PT::Array{Array{Array{Int, 1}, 1}, 1}
#	- output2 from ysr()
# Return Values:
#	FFT::Array{BigInt, 2}
#	- FFT is the Fast Fourier Transform of SNF
function sn_fft(N::Int, M::Int, SNF::Array{BigInt, 1}, YSR::Array{Array{Array{SparseMatrixCSC, 1}, 1}, 1}, PT::Array{Array{Array{Int, 1}, 1}, 1})
	np = nprocs()
	if np == 1  # || N < 6
		C = Counter(1)
		return compute_fft(N, M, SNF, YSR, PT, C)
	else #This is a modified version of Julia's pmap()
		sFFT = Array(Array{Array{BigInt, 2}, 1}, N)
		RR_YSR = Array(RemoteRef, np)
		RR_PT = Array(RemoteRef, np)
		for p = 1:np
			if p != myid()
				RR_YSR[p] = RemoteRef(p)
				put!(RR_YSR[p], YSR)
				RR_PT[p] = RemoteRef(p)
				put!(RR_PT[p], PT)
			end
		end
		BS = factorial(N - 1)
		i = 1
		nextidx() = (idx = i; i += 1; idx)
		@sync begin
			for p = 1:np
				if p != myid()
					@async begin
						while true
							idx = nextidx()
							if idx > N
								break
							end
							lb = (idx - 1) * BS + 1
							ub = lb + BS - 1
							SNF_subgroup = SNF[lb:ub]
							sFFT[idx] = remotecall_fetch(p, compute_fft_remote, N - 1, SNF_subgroup, RR_YSR[p], RR_PT[p], Counter(1))
						end
					end
				end
			end
		end
		YSRn = YSR[N]
		NP = length(YSRn)
		PTn = PT[N]
		RR_sFFT = Array(RemoteRef, np)
		for p = 1:np
			if p != myid()
				RR_sFFT[p] = RemoteRef(p)
				put!(RR_sFFT[p], sFFT)
			end
		end
		FFT = Array(Array{BigInt, 2}, NP)
		i = 1
		nextidx() = (idx = i; i += 1; idx)
		@sync begin
			for p = 1:np
				if p != myid()
					@async begin
						while true
							idx = nextidx()
							if idx > NP
								break
							end
							FFT[idx] = remotecall_fetch(p, fc_remote, N, idx, RR_YSR[p], RR_PT[p], RR_sFFT[p])
						end
					end
				end
			end
		end
		return FFT
	end
end

# Parameters:
#	N::Int
#	- the size of the FFT being calculated is N
#       M::Int
#       - the module
#	SNF::Array{BigInt, 1}
#	- SNF[i] is the value associated with the Permutation that permutation_index() maps to i
#	YSR::Array{Array{Array{SparseMatrixCSC, 1}, 1}, 1}
#	- output1() from ysr()
#	PT::Array{Array{Array{Int, 1}, 1}, 1}
#	- output2() from ysr()
#	C::Counter
#	- a wrapper for an Int that is used to coordinate the recursive access to the elements of SNF
# Return Values:
#       FFT::Array{Array{BigInt, 2}, 1}
#	- FFT is the Fast Fourier Transfrom of SNF
function compute_fft(N::Int, M::Int, SNF::Array{BigInt, 1}, YSR::Array{Array{Array{SparseMatrixCSC, 1}, 1}, 1}, PT::Array{Array{Array{Int, 1}, 1}, 1}, C::Counter)
	if N == 1
		sFFT = Array(Array{BigInt, 2}, 1)
		sFFTi = Array(BigInt, 1, 1)
		sFFTi[1, 1] = SNF[C.N]
		sFFT[1] = sFFTi
		C.N += 1
		return sFFT
	end
	sFFT = Array(Array{Array{BigInt, 2}, 1}, N)
	for n = 1:N
		sFFT[n] = compute_fft(N - 1, M, SNF, YSR, PT, C)
	end
	FFT = combine_sfft(N, M, YSR[N], PT[N], sFFT)
	return FFT
end

function compute_fft_remote(N::Int, M::Int, SNF::Array{BigInt, 1}, RR_YSR::RemoteRef, RR_PT::RemoteRef, C::Counter)
	YSR = fetch(RR_YSR)
	PT = fetch(RR_PT)
	FFT =  compute_fft(N, M, SNF, YSR, PT, C)
	return FFT
end

# Parameters:
#	N::Int
#	- the size of the FFT being calculated is N
#       M::Int
#       - the module
#	YSRn::Array{Array{SparseMatrixCSC, 1}, 1}, 1}
#	- Young's Seminormal Representations for the Partitions of N
#	PTn::Array{Array{Int, 1}, 1}
#	- the decomposition indices for the Partitions of N
#	sFFT::Array{Array{Array{BigInt, 2}, 1}, 1}
#	- the array of N FFT's of size N-1 that is used to compute this FFT of size N
# Return Values:
#	FFT::Array{Array{BigInt, 2}, 1}
#	- FFT is the Fast Fourier Transform for this group
function combine_sfft(N::Int, M::Int, YSRn::Array{Array{SparseMatrixCSC, 1}, 1}, PTn::Array{Array{Int, 1}, 1}, sFFT::Array{Array{Array{BigInt, 2}, 1}, 1})
	NP = length(YSRn)
	FFT = Array(Array{BigInt, 2}, NP)
	for p = 1:NP
		YSRnp = YSRn[p]
		PTnp = PTn[p]
		FC = fc(N, M, YSRnp, PTnp, sFFT)
		FFT[p] = FC
	end
	return FFT
end

# Parameters:
#	N::Int
#	- the Fourier Component that is being calculated corresponds to a Partition, P, of N
#       M::Int
#       - the module
#	YSRnp::Array{SparseMatrixCSC, 1}
#	- Young's Orthogonal Representations for P
#	PTnp::Array{Int, 1}
#	- the indices for the Partitions that P decomposes into
#	sFFT::Array{Array{Array{BigInt, 2}, 1}, 1}
#	- the array of N FFT's of size N-1 that is used to compute this FFT of size N
# Return Values:
#	FC::Array{BigInt, 2}
#	- FC is the Fourier Coefficient corresponding to the Partition P
function fc(N::Int, M::Int, YSRnp::Array{SparseMatrixCSC, 1}, PTnp::Array{Int, 1}, sFFT::Array{Array{Array{BigInt, 2}, 1}, 1})
	Dim = size(YSRnp[1], 1)
	FC = dsm(Dim, sFFT[N], PTnp)
	CCM = eye(BigInt,Dim)
	for n = (N - 1):-1:1
    # compute representation associated with [[n,N]]
		CCM = (big(YSRnp[n]) * CCM) % M
		DSM = big(dsm(Dim, sFFT[n], PTnp)) % M
		FC_n = (big(CCM) * big(DSM)) % M
		FC = (big(FC) + FC_n) % M
	end
	return (FC % M + M) % M
end

function fc_remote(N::Int, M::Int, p::Int, RR_YSR::RemoteRef, RR_PT::RemoteRef, RR_sFFT::RemoteRef)
	YSR = fetch(RR_YSR)
	PT = fetch(RR_PT)
	sFFT = 	fetch(RR_sFFT)
	FC = fc(N, M, YSR[N][p], PT[N][p], sFFT)
	return FC
end

# Parameters:
#	Dim::Int
#	- the size of the DSM that is being caculated
#	sFFTn::Array{Array{BigInt, 2}, 1}
#	- one of the elements of sFFT from fc()
#	PTnp::Array{Int, 1}
#	- same as in fc()
# Return Values:
#	DSM::Array{BigInt, 2} (Direct Sum Matrix)
#	- the direct sum of the coefficients of sFFTn that correspond to the Partitions indicated by PTnp
function dsm(Dim::Int, sFFTn::Array{Array{BigInt, 2}, 1}, PTnp::Array{Int, 1})
	DSM = zeros(BigInt, Dim, Dim)
	offset = 0
	for i = 1:length(PTnp)
		FC = sFFTn[PTnp[i]]
		sDim = size(FC, 1)
		for r = 1:sDim
			for c = 1:sDim
				v = FC[r, c]
				if v != 0
					DSM[offset + r, offset + c] = v
				end
			end
		end
		offset += sDim
	end
	return DSM
end
