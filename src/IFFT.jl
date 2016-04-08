# Parameters:
#	N::Int
#	- the problem size
#	M::Int
#	- the module
#	FFT::Array{Array{BigInt, 2}, 1}
#	- a Fast Fourier Transform of size N
#	- should be the output of sn_fft() or sn_fft_sp()
#	YSR::Array{Array{Array{SparseMatrixCSC, 1}, 1}, 1}
#	- output1 from YSR()
#	PT::Array{Array{Array{Int, 1}, 1}, 1}
#	- output2 from ysr()
# Return Values:
#	SNF::Array{BigInt, 1}
#	- the function over Sn that corresponds to FFT
function sn_ifft(N::Int, M::Int, FFT::Array{Array{BigInt, 2}, 1},  YSR::Array{Array{Array{SparseMatrixCSC, 1}, 1}, 1}, PT::Array{Array{Array{Int, 1}, 1}, 1})
	np = nprocs()
	if np == 1 || N < 6
		SNF = Array(BigInt, factorial(N))
		compute_ifft(N, M, SNF, FFT, YSR, PT, Counter(1))
		return SNF
	else
		pSNFA = Array(Array{BigInt, 1}, N)
		RR_FFT = Array(RemoteRef, np)
		RR_YSR = Array(RemoteRef, np)
		RR_PT = Array(RemoteRef, np)
		for p = 1:np
			if p != myid()
				RR_FFT[p] = RemoteRef(p)
				put!(RR_FFT[p], FFT)
				RR_YSR[p] = RemoteRef(p)
				put!(RR_YSR[p], YSR)
				RR_PT[p] = RemoteRef(p)
				put!(RR_PT[p], PT)
			end
		end
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
							pSNFA[idx] = remotecall_fetch(p, compute_sifft_remote, N, idx, RR_FFT[p], RR_YSR[p], RR_PT[p])
						end
					end
				end
			end
		end
		SNF = Array(BigInt, factorial(N))
		BS = factorial(N - 1)
		i = 1
		for n = 1:N
			pSNF = pSNFA[n]
			for si = 1:BS
				SNF[i] = pSNF[si]
				i += 1
			end
		end
		return SNF
	end
end

# Parameters:
#	N::Int
#	- the problem size
#	M::Int
#	- the module
#	SNF::Array{BigInt, 1}
#	- the function over Sn that is being calculated
#	FFT::Array{Array{BigInt, 2}, 1}
#	- a Fast Fourier Transform of size N
#	- should be the output of sn_fft() or sn_fft_sp()
#	YSR::Array{Array{Array{SparseMatrixCSC, 1}, 1}, 1}
#	- output1 from ysr()
#	PT::Array{Array{Array{Int, 1}, 1}, 1}
#	- output2 from ysr()
#	C::Counter
#	- a wrapper for an Int that is used to coordinate the recursive defining of the elements of SNF
# Return Values:
#	None
#	- once the recursion is finished, SNF is full
function compute_ifft(N::Int, M::Int, SNF::Array{BigInt, 1}, FFT::Array{Array{BigInt, 2}, 1}, YSR::Array{Array{Array{SparseMatrixCSC, 1}, 1}, 1}, PT::Array{Array{Array{Int, 1}, 1}, 1}, C::Counter)
	if N == 2
		YSRn = YSR[2]
		sFFT = (invmod(big(2),M) * ((big(full(YSRn[1][1])) * FFT[1]) %  M + (big(full(YSRn[2][1])) * FFT[2])) % M) % M
		SNF[C.N] = sFFT[1,1]
		C.N += 1
		sFFT = (invmod(big(2),M) * (FFT[1] + FFT[2]) % M) % M
		SNF[C.N] = sFFT[1,1]
		C.N += 1
	else
		YSRn = YSR[N]
		NPn = length(YSRn)
		YSRd = YSR[N - 1]
		NPd = length(YSRd)
		PTn = PT[N]
		sFFT = Array(Array{BigInt, 2}, NPd)
		for n = 1:N
			for p = 1:NPd
				Dim = size(YSRd[p][1], 1)
				sFFT[p] = zeros(BigInt, Dim, Dim)
			end
			for p = 1:NPn
				update_sfft!(N, M, n, sFFT, FFT[p], YSRn[p], YSRd, PTn[p])
			end
			compute_ifft(N - 1, M, SNF, sFFT, YSR, PT, C)
		end
	end
end

# Parameters:
#	N::Int
#	- the problem size
#	M::Int
#	- the module
#	n::Int
#	- determines which left-sided coset we are defining the FFT for
#	sFFT::Array{Array{BigInt, 2}, 1}
#	- the FFT of the nth left-sided coset that we are defining
#	FFTp::Array{BigInt, 2}
#	- the pth component of the FFT of size N
#	YSRnp::Array{SparseMatrixCSC, 1}, 1}
#	- Youngs Seminormal Representations for the pth Partition of N
#	YSRd::Array{Array{SparseMatrixCSC, 1}, 1}
#	- Youngs Seminormal Representations for the Partitions of N - 1
#	PTnp::Array{Array{Array{Int, 1}, 1}, 1}
#	- The decomposition indices for the pth Partition of N
# Return Values:
#	None
#	- sFFT is fully defined once all of the components of FFT have been used
function update_sfft!(N::Int, M::Int, n::Int, sFFT::Array{Array{BigInt, 2}, 1}, FFTp::Array{BigInt, 2}, YSRnp::Array{SparseMatrixCSC, 1}, YSRd::Array{Array{SparseMatrixCSC, 1}, 1}, PTnp::Array{Int, 1})
        Dim = size(YSRnp[1], 1)
	P = eye(BigInt, Dim, Dim)
  # compute the inverse of [[n,N]]
	for ccn = n:(N - 1)
	  P = (big(YSRnp[ccn]) * P) % M # [original SnFFT] P = P * YSRnp[ccn]
	end
  # [original SnFFT] P = transpose(P)
  P = P % M
	P = (P * big(FFTp)) % M
	lb = 1
	for d = 1:length(PTnp)
		index = PTnp[d]
		sDim = size(YSRd[index][1], 1)
		ub = lb + sDim - 1
		sFFT[index] = (sFFT[index] + (M + (Dim * invmod(big(sDim * N), M)) * big(P[lb:ub, lb:ub]))) % M
		lb += sDim
	end
end

# Parameters:
#	N::Int
#	- the problem size
#	M::Int
#	- the module
#	n::Int
#	- determines which left-sided coset we are defining the function for
#	FFT::Array{Array{BInt, 2}, 1}
#	- a Fast Fourier BigTransform of size N
#	- should be the output of sn_fft() or sn_fft_sp()
#	YSR::Array{Array{Array{SparseMatrixCSC, 1}, 1}, 1}
#	- output1 from ysr()
#	PT::Array{Array{Array{Int, 1}, 1}, 1}
#	- output2 from ysr()
# Return Values:
#	pSNF::Array{BigInt, 1}
#	- the function on the nth left-sided coset of Sn
function compute_sifft(N::Int, M::Int, n::Int, FFT::Array{Array{Int, 2}, 1}, YSR::Array{Array{Array{SparseMatrixCSC, 1}, 1}, 1},  PT::Array{Array{Array{Int, 1}, 1}, 1})
	pSNF = Array(BigInt, factorial(N - 1))
	C = Counter(1)
	YSRn = YSR[N]
	NPn = length(YSRn)
	YSRd = YSR[N - 1]
	NPd = length(YSRd)
	PTn = PT[N]
	sFFT = Array(Array{BigInt, 2}, NPd)
	for p = 1:NPd
		Dim = size(YSRd[p][1], 1)
		sFFT[p] = zeros(BigInt, Dim, Dim)
	end
	for p = 1:NPn
		update_sfft!(N, M, n, sFFT, FFT[p], YSRn[p], YSRd, PTn[p])
	end
	compute_ifft(N - 1, M, pSNF, sFFT, YSR, PT, C)
	return pSNF
end

function compute_sifft_remote(N::Int, M::Int, n::Int, RR_FFT::RemoteRef, RR_YSR::RemoteRef, RR_PT::RemoteRef)
	FFT = fetch(RR_FFT)
	YSR = fetch(RR_YSR)
	PT = fetch(RR_PT)
	pSNF = compute_sifft(N, M, n, FFT, YSR, PT)
	return pSNF
end
