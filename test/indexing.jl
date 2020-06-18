using TriMatrices: trinum, triinv, triinv_rem
using TriMatrices: mat2tri, tri2mat


@testset "tri and triinv" begin
    for n in 0:100
        t = trinum(n)

        @test t == sum(1:n)

        @test triinv(t) == n
        @test triinv(t, true) == n
    end

    for t in 0:100
        n, r = triinv_rem(t)
        @test n >= 0 && r >= 0
        @test tri(n) + r == t
        @test triinv(t) == n
        r == 0 || @test_throws DomainError triinv(t, true)
    end
end


@testset "mat2tri and tri2mat" begin
    i = 1
    for row in 1:5
        for col in 1:row
            @test mat2tri(row, col) == i
            @test tri2mat(i) == (row, col)
            i += 1
        end
    end

    @test_throws DomainError mat2tri(1, 2)
    @test_throws DomainError mat2tri(0, 0)
end
