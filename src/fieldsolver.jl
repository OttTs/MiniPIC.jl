#=
Discretization nodes
          1   2   3   4   5
        φ |---|---|---|---|---
MC-PIC  E |---|---|---|---|---
EC-PIC  E --|---|---|---|---|–
            1   2   3   4   5

After non-dimensionalization, Δx = 1
=#

function fieldsolver!(x, E, α, do_MC_PIC::Bool)
    Nₚ = length(x)
    Nₓ = length(E)

    # Deposit charge to the grid ρ = ∑S(x - xₚ)
    E .= - Nₚ / Nₓ  # Neutralizing background
    for xₚ in x
        deposit!(xₚ, E)
    end

    # Solve Poisson equation and get electric field
    get_potential!(E)
    get_field!(E, do_MC_PIC)
    E .*= α
end


function get_potential!(φ)
    # Poisson solver φ = Δ⁻¹ ρ (Thomas algorithm)
    # φ must be filled with ρ
    Nₓ = length(φ)

    # Forward elimination
    φ[1] = 0
    for iₓ in 2:Nₓ
        φ[iₓ] = φ[iₓ] + (iₓ - 1) / iₓ * φ[iₓ-1]

    end

    # Back substitution
    φ[end] = - Nₓ /(Nₓ + 1) * φ[end]
    for iₓ in (Nₓ - 1):-1:1
        φ[iₓ] = iₓ / (iₓ + 1) * (φ[iₓ+1] - φ[iₓ])
    end
end


function get_field!(E, do_MC_PIC::Bool)
    # TODO Theoretically, this should be -E (Due to nondimensionalization)
    # E must be filled with φ
    Nₓ = length(E)

    # Calculate the staggered electric field
    # Eᵢ = φᵢ₊₁ - φᵢ / Δxᵢ
    φ₁ = E[1]
    for iₓ in 1:(Nₓ - 1)
        E[iₓ] = E[iₓ+1] - E[iₓ]
    end
    E[Nₓ] = φ₁ - E[Nₓ]

    if do_MC_PIC
        # For MC-PIC, convert to centered nodes
        # Eᵢ = (Eᵢ + Eᵢ₋₁) / 2 = (φᵢ₊₁ - φᵢ₋₁) / 2Δx
        Eₙ = E[Nₓ]
        for iₓ in (Nₓ:-1:2)
            E[iₓ] = 0.5(E[iₓ] + E[iₓ-1])
        end
        E[1] = 0.5(E[1] + Eₙ)
    end
end


function interpolate(xₚ, E, do_MC_PIC::Bool)
    # Interpolate the electric field to particle positions
    Nₓ = length(E)
    iₓ = cellindex(xₚ, Nₓ)

    if do_MC_PIC
        # MC-PIC: linear interpolation at centered nodes
        iₓ₊₁ = mod(iₓ, Nₓ) + 1
        ω = xₚ - iₓ + 1
        return (1 - ω) * E[iₓ] + ω * E[iₓ₊₁]
    else
        # EC-PIC: constant interpolation at staggered nodes
        return E[iₓ]
    end
end

function deposit!(xₚ, ρ)
    # Deposit particle charge to the grid
    Nₓ = length(ρ)
    iₓ = cellindex(xₚ, Nₓ)
    iₓ₊₁ = mod(iₓ, Nₓ) + 1
    ω = xₚ - iₓ + 1
    ρ[iₓ] += (1 - ω)
    ρ[iₓ₊₁] += ω
end

cellindex(xₚ, Nₓ) = max(ceil(Int, xₚ), 1)