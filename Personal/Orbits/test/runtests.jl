# test/runtests.jl

using Orbits
using Test
using OrdinaryDiffEq


@testset "Orbits.jl Tests" begin

    @testset "Body Constructor" begin
        # Default constructor
        b1 = Body()
        @test b1.name == "Body"
        @test b1.mass == 1.0
        @test b1.x == 0.0
        @test b1.y == 0.0
        @test b1.vx == 0.0
        @test b1.vy == 0.0
        @test b1.color == :black
        @test b1.size == 15.0
        
        # Custom constructor
        b2 = Body(name="Earth", mass=3e-6, pos=(1.0, 0.0), vel=(0.0, 1.0), color=:blue, size=8.0)
        @test b2.name == "Earth"
        @test b2.mass == 3e-6
        @test b2.x == 1.0
        @test b2.y == 0.0
        @test b2.vx == 0.0
        @test b2.vy == 1.0
        @test b2.color == :blue
        @test b2.size == 8.0
    end

    @testset "Kepler Orbit Physics & Conservation" begin
        # Set up a circular orbit of a light body (Earth, m = 1e-6) around a heavy body (Sun, M = 1.0)
        # R = 1.0, G = 1.0. Circular orbit velocity: v = sqrt(G*M/R) = 1.0.
        # Theoretical orbital period: T = 2*pi*R / v = 2*pi.
        
        sun = Body(name="Sun", mass=1.0, pos=(0.0, 0.0), vel=(0.0, 0.0), color=:yellow, size=30)
        earth = Body(name="Earth", mass=1e-6, pos=(1.0, 0.0), vel=(0.0, 1.0), color=:blue, size=10)
        
        bodies = [sun, earth]
        G = 1.0
        T = 2 * pi
        
        # Solve with high-accuracy solver and tight tolerances
        sol = simulate_orbits(bodies, (0.0, T); G=G, solver=Tsit5(), reltol=1e-10, abstol=1e-10)

        
        # Test 1: Orbits evolved to the specified end time
        @test sol.t[end] ≈ T
        
        # Helper function to compute total energy
        function compute_energy(state, masses, G)
            N = length(masses)
            K = 0.0 # Kinetic energy
            V = 0.0 # Potential energy
            
            # Kinetic energy
            for i in 1:N
                vx = state[4*(i-1) + 3]
                vy = state[4*(i-1) + 4]
                K += 0.5 * masses[i] * (vx^2 + vy^2)
            end
            
            # Gravitational potential energy
            for i in 1:N
                xi = state[4*(i-1) + 1]
                yi = state[4*(i-1) + 2]
                for j in (i+1):N
                    xj = state[4*(j-1) + 1]
                    yj = state[4*(j-1) + 2]
                    r = sqrt((xj - xi)^2 + (yj - yi)^2)
                    V -= G * masses[i] * masses[j] / r
                end
            end
            
            return K + V
        end
        
        masses = [b.mass for b in bodies]
        
        E_init = compute_energy(sol.u[1], masses, G)
        E_final = compute_energy(sol.u[end], masses, G)
        
        # Test 2: Energy conservation (within very tight numerical tolerance)
        @test E_final ≈ E_init atol=1e-9
        
        # Check energy at intermediate points
        for t in range(0.0, T, length=10)
            E_t = compute_energy(sol(t), masses, G)
            @test E_t ≈ E_init atol=1e-9
        end
        
        # Test 3: Return to initial conditions after 1 period
        state_final = sol.u[end]
        
        # Sun position and velocity should remain close to 0
        @test state_final[1] ≈ 0.0 atol=1e-4 # Sun x
        @test state_final[2] ≈ 0.0 atol=1e-4 # Sun y
        @test state_final[3] ≈ 0.0 atol=1e-4 # Sun vx
        @test state_final[4] ≈ 0.0 atol=1e-4 # Sun vy
        
        # Earth position and velocity should return to (1,0) and (0,1)
        @test state_final[5] ≈ 1.0 atol=1e-4 # Earth x
        @test state_final[6] ≈ 0.0 atol=1e-4 # Earth y
        @test state_final[7] ≈ 0.0 atol=1e-4 # Earth vx
        @test state_final[8] ≈ 1.0 atol=1e-4 # Earth vy
    end

end
