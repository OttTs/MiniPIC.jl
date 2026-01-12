using MiniPIC
using Test

@testset "MiniPIC.jl" begin
    @testset "Field Solver Internals" begin
        # Coordinate system: domain is [0, Nₓ) where Nₓ is number of cells
        # Cell i spans (i-1, i] for i = 1, 2, ..., Nₓ (using ceil)

        @testset "cellindex" begin
            Nₓ = 10
            # cellindex uses ceil, so cell 1 covers [0, 1], cell 2 covers (1, 2], etc.
            @test MiniPIC.cellindex(0.0, Nₓ) == 1   # At origin, maps to cell 1
            @test MiniPIC.cellindex(0.5, Nₓ) == 1   # In first cell
            @test MiniPIC.cellindex(1.0, Nₓ) == 1   # At boundary (ceil(1)=1, max(1,1)=1)
            @test MiniPIC.cellindex(1.1, Nₓ) == 2   # Just past boundary
            @test MiniPIC.cellindex(5.5, Nₓ) == 6   # Middle of domain
            @test MiniPIC.cellindex(9.9, Nₓ) == 10  # Near end
        end

        @testset "deposit!" begin
            # Test charge deposition with linear weighting
            # Domain: [0, 4), particles deposit to nodes 1-4
            Nₓ = 4

            # Particle at x=1.5 (cell 2, weight splits between nodes 2 and 3)
            ρ = zeros(Nₓ)
            MiniPIC.deposit!(1.5, ρ)
            @test ρ[2] ≈ 0.5  # Half weight to node 2
            @test ρ[3] ≈ 0.5  # Half weight to node 3
            @test ρ[1] ≈ 0.0
            @test ρ[4] ≈ 0.0

            # Particle at node boundary (x=2.0)
            ρ = zeros(Nₓ)
            MiniPIC.deposit!(2.0, ρ)
            @test ρ[3] ≈ 1.0  # Full weight at node 3
            @test ρ[1] ≈ 0.0
            @test ρ[2] ≈ 0.0

            # Particle in cell 4 (near periodic boundary, x=3.5)
            ρ = zeros(Nₓ)
            MiniPIC.deposit!(3.5, ρ)
            @test ρ[4] ≈ 0.5
            @test ρ[1] ≈ 0.5  # Wraps to node 1 (periodic)
        end

        @testset "interpolate MC-PIC" begin
            # Test linear interpolation at centered nodes
            # For MC-PIC with E = [1, 2, 3, 4]:
            E = [1.0, 2.0, 3.0, 4.0]

            # At x=0.0: cell 1, ω = 0-1+1 = 0, result = (1-0)*E[1] + 0*E[2] = 1.0
            @test MiniPIC.interpolate(0.0, E, true) ≈ 1.0

            # At x=0.5: cell 1, ω = 0.5, result = 0.5*E[1] + 0.5*E[2] = 1.5
            @test MiniPIC.interpolate(0.5, E, true) ≈ 1.5

            # At x=1.0: cell 1, ω = 1.0, result = 0*E[1] + 1*E[2] = 2.0
            @test MiniPIC.interpolate(1.0, E, true) ≈ 2.0

            # At x=1.5: cell 2, ω = 0.5, result = 0.5*E[2] + 0.5*E[3] = 2.5
            @test MiniPIC.interpolate(1.5, E, true) ≈ 2.5
        end

        @testset "interpolate EC-PIC" begin
            # Test constant interpolation at staggered nodes
            E = [1.0, 2.0, 3.0, 4.0]

            # EC-PIC uses constant interpolation within each cell
            @test MiniPIC.interpolate(0.0, E, false) ≈ 1.0  # Cell 1
            @test MiniPIC.interpolate(0.5, E, false) ≈ 1.0  # Cell 1
            @test MiniPIC.interpolate(1.5, E, false) ≈ 2.0  # Cell 2
            @test MiniPIC.interpolate(2.5, E, false) ≈ 3.0  # Cell 3
        end

        @testset "get_potential! (Poisson solver)" begin
            # Test that Poisson solver produces finite results
            Nₓ = 10
            ρ = ones(Nₓ)
            MiniPIC.get_potential!(ρ)
            @test all(isfinite, ρ)

            # Test with sinusoidal charge density
            ρ_sin = [sin(2π * i / Nₓ) for i in 1:Nₓ]
            MiniPIC.get_potential!(ρ_sin)
            @test all(isfinite, ρ_sin)
        end

        @testset "get_field!" begin
            # Test that field computation produces finite results
            Nₓ = 4
            φ = collect(1.0:Nₓ)

            # MC-PIC (centered nodes)
            E_mc = copy(φ)
            MiniPIC.get_field!(E_mc, true)
            @test all(isfinite, E_mc)

            # EC-PIC (staggered nodes)
            E_ec = copy(φ)
            MiniPIC.get_field!(E_ec, false)
            @test all(isfinite, E_ec)

            # For constant potential, field should be zero
            φ_const = ones(Nₓ)
            E_zero = copy(φ_const)
            MiniPIC.get_field!(E_zero, true)
            @test all(E -> abs(E) < 1e-10, E_zero)
        end
    end

    @testset "Energy Calculations" begin
        @testset "kinetic_energy" begin
            # Single particle at rest
            @test MiniPIC.kinetic_energy([0.0], 1.0) ≈ 0.0

            # Single particle with velocity 1: E = 0.5 * m * v^2 = 0.5
            @test MiniPIC.kinetic_energy([1.0], 1.0) ≈ 0.5

            # Two particles: E = 0.5 * m * (v1^2 + v2^2) = 0.5 * (1 + 4) = 2.5
            @test MiniPIC.kinetic_energy([1.0, 2.0], 1.0) ≈ 2.5

            # With different mass: E = 0.5 * 3 * 4 = 6
            @test MiniPIC.kinetic_energy([2.0], 3.0) ≈ 6.0
        end

        @testset "potential_energy" begin
            # Zero field
            @test MiniPIC.potential_energy([0.0, 0.0], 1.0, 1.0) ≈ 0.0

            # Uniform field: E_pot = 0.5 * ε * sum(E^2) * Δx
            # With E = [1,1,1], Δx = 0.5, ε = 2: E_pot = 0.5 * 2 * 3 * 0.5 = 1.5
            @test MiniPIC.potential_energy([1.0, 1.0, 1.0], 0.5, 2.0) ≈ 1.5
        end
    end

    @testset "Fieldsolver Integration" begin
        @testset "Uniform distribution" begin
            # Uniformly distributed particles should produce near-zero net field
            Nₚ = 1000
            Nₓ = 10
            x = collect(range(0.0, Nₓ - Nₓ/Nₚ, length=Nₚ))
            E = zeros(Nₓ)
            α = 1.0

            MiniPIC.fieldsolver!(x, E, α, true)

            # For uniform distribution, electric field should be small
            @test maximum(abs, E) < 0.5
        end

        @testset "Charge conservation" begin
            # Total deposited charge should equal number of particles
            Nₚ = 100
            Nₓ = 10
            x = rand(Nₚ) .* Nₓ
            ρ = zeros(Nₓ)

            for xₚ in x
                MiniPIC.deposit!(xₚ, ρ)
            end

            @test sum(ρ) ≈ Nₚ
        end

        @testset "MC-PIC vs EC-PIC" begin
            # Both methods should produce finite results
            Nₚ = 100
            Nₓ = 10
            x = rand(Nₚ) .* Nₓ
            α = 1.0

            E_mc = zeros(Nₓ)
            E_ec = zeros(Nₓ)

            MiniPIC.fieldsolver!(copy(x), E_mc, α, true)
            MiniPIC.fieldsolver!(copy(x), E_ec, α, false)

            @test all(isfinite, E_mc)
            @test all(isfinite, E_ec)
        end
    end

    @testset "Simulation Steps" begin
        @testset "Free streaming" begin
            # With α=0 (no coupling), particles should move freely
            Nₓ = 10
            x = [5.0]
            v = [0.1]
            E = zeros(Nₓ)
            α = 0.0

            v_init = v[1]
            MiniPIC.simulate_steps!(x, v, E, α; num_steps=10, do_MC_PIC=true, progress=nothing)

            # Velocity change should be minimal
            @test abs(v[1] - v_init) < 1.0
        end

        @testset "Periodic boundary conditions" begin
            # Particles should wrap around at domain boundaries
            Nₓ = 10
            x = [9.5]   # Near right boundary
            v = [1.0]   # Moving right
            E = zeros(Nₓ)
            α = 0.0

            MiniPIC.simulate_steps!(x, v, E, α; num_steps=1, do_MC_PIC=true, progress=nothing)

            # Position should be in valid range after wrapping
            @test 0.0 <= x[1] < Nₓ
        end

        @testset "Multiple particles" begin
            # Test with multiple particles
            Nₓ = 10
            Nₚ = 50
            x = rand(Nₚ) .* Nₓ
            v = randn(Nₚ) .* 0.1
            E = zeros(Nₓ)
            α = 0.1

            # Run some steps
            MiniPIC.simulate_steps!(x, v, E, α; num_steps=5, do_MC_PIC=true, progress=nothing)

            # All positions should be valid
            @test all(0.0 .<= x .< Nₓ)
            @test all(isfinite, v)
        end

        @testset "Reversibility check" begin
            # Leapfrog is symplectic - test numerical stability
            Nₓ = 10
            Nₚ = 20
            x = rand(Nₚ) .* Nₓ
            v = randn(Nₚ) .* 0.1
            E = zeros(Nₓ)
            α = 0.05

            # Run forward
            MiniPIC.simulate_steps!(x, v, E, α; num_steps=10, do_MC_PIC=true, progress=nothing)

            # Check results are finite
            @test all(isfinite, x)
            @test all(isfinite, v)
        end
    end

    @testset "Physical Scenarios" begin
        @testset "Two-stream setup" begin
            # Two counter-streaming beams should have zero net momentum
            Nₚ = 100
            v_beam = 1.0

            x = rand(Nₚ) .* 10
            v = vcat(fill(v_beam, Nₚ÷2), fill(-v_beam, Nₚ÷2))

            @test length(x) == Nₚ
            @test length(v) == Nₚ
            @test sum(v)/length(v) ≈ 0.0 atol=0.01
        end

        @testset "Energy behavior" begin
            # Total energy should remain bounded during simulation
            Nₓ = 20
            Nₚ = 100
            x = collect(range(0.0, Nₓ - Nₓ/Nₚ, length=Nₚ))
            v = randn(Nₚ) .* 0.1
            E = zeros(Nₓ)
            α = 0.01

            E_kin_init = MiniPIC.kinetic_energy(v, 1.0)

            MiniPIC.simulate_steps!(x, v, E, α; num_steps=5, do_MC_PIC=true, progress=nothing)

            E_kin_final = MiniPIC.kinetic_energy(v, 1.0)

            # Energy should not blow up
            @test E_kin_final < 10 * E_kin_init
            @test E_kin_final > 0
        end

        @testset "Cold plasma" begin
            # Uniform cold plasma should remain stable
            Nₓ = 10
            Nₚ = 100
            x = collect(range(0.0, Nₓ - Nₓ/Nₚ, length=Nₚ))
            v = zeros(Nₚ)  # Cold plasma (no thermal velocity)
            E = zeros(Nₓ)
            α = 0.1

            MiniPIC.simulate_steps!(x, v, E, α; num_steps=10, do_MC_PIC=true, progress=nothing)

            # Velocities should remain small
            @test maximum(abs, v) < 1.0
        end
    end

end
