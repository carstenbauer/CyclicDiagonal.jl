using CyclicDiagonal
using Test
using SpecialMatrices
using LinearAlgebra

@testset "CyclicDiagonal.jl" begin
    # Write your own tests here.

    M = Matrix( Circulant([1;3;2]) )

    # Test simple 3x3 circulant matric
    D1, C, D2 = CyclicDiagonal.iterDCD( M )
    @test norm(D1-I, Inf) < 1e-12

    # Test simple product of a diagonal and circulant matrix
    D = Diagonal(1:3)
    D1, C, D2 = CyclicDiagonal.iterDCD( D*M )
    @test norm( D1*C*D2 - D*M, Inf) < 1e-12

    # Test products of kxk matrices (k=1,...,K), each size N times
    N = 10
    K = 10
    r = 1:50
    for k in 1:K
        for i in 1:N
            D1o = Diagonal( rand(r, k) )
            Co  = Circulant( rand(r, k) )
            D2o = Diagonal( rand(r, k) )
            D1, C, D2 = CyclicDiagonal.iterDCD( D1o * Co * D2o, 1e6, 1e-7 )
            @test norm( D1*C*D2 - D1o*Co*D2o, Inf) < 1e-3
        end
    end
end
