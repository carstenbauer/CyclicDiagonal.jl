module CyclicDiagonal

using SpecialMatrices, LinearAlgebra

export iterDCD

# Iteratively find an approximate factorization of A as a product D1 * C * D2
# where D1 and D2 are diagonal matrices, and C is a circulant matrix
function iterDCD(A::Matrix, max_iter = 1e5, tol = 1e-5)
	s_col = ones(size(A,1))
	s_row = Array{Float64}(undef, 0, 2)
	mu_row = Array{Float64}(undef, 0, 2)
	mu_col = Array{Float64}(undef, 0, 2)

	d = Inf
	i = 0
	while i < max_iter && d > tol
		mu_row, s_row = find_mu_and_s( A .* (s_col.^-1)', max_iter, tol)
		mu_col, s_col = find_mu_and_s( ((s_row.^-1) .* A)', max_iter, tol )
		d = norm( [ mu_col[1]; mu_col[end:-1:2] ]  - mu_row )
		i = i + 1
	end

	D1 = Diagonal( s_row[:] )
	C =  Circulant( mu_col[:] )
	D2 = Diagonal( s_col[:] )

	D1, C, D2
end

# Find left side diagonal scaling matrix with diagonal s and circulant diagonal matrix factorization approximations
function find_mu_and_s(A, max_iter , tol)

    uncirc!(A)
    mu = sum(A, dims=1)'/size(A,1)
    mu_old = fill(Inf, size(mu))
	s = ones(size(A,1))

    i = 0
    while i < max_iter && norm( mu_old - mu ) > tol

        s = A*mu/(mu'*mu)
        mu_old = mu
        mu = sum(A ./ s, dims=1)'/size(A,1)

        i = i + 1
    end

	mu, s
end

# Circulate rows of A
function circ!(A)
    for i in 1:size(A,1)
        A[i,:] = circshift(A[i,:], (i-1))
    end
	A
end

# Uncirculate rows of A
function uncirc!(A)
    for i in 1:size(A,1)
        A[i,:] = circshift(A[i,:], (1-i))
    end
	A
end

end # module
