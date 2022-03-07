using LinearAlgebra
using StaticArrays
using BenchmarkTools
using FiniteDiff
using Test
using EndoBeams: Sphere_SDF, Plane_z_SDF, contact_gap, local_R⁰, rotation_matrix, compute, ID3, Tₛ⁻¹







# # information from node 1 and 2
const X₁ = @SVector [0.91, 0.58, 0.68]
const X₂ = @SVector [0.97, 0.04, 0.41]

const l₀ = norm(X₂-X₁)

u₁ = @SVector [0.8, 0.5, 0.1]
u₂ = @SVector [0.9, 0.52, 0.14]

Θ₁ = @SVector [-0.2423, 0.9047, 0.5896]
Θ₂ = @SVector [-0.2394, 0.8994, 0.5902]


u̇₁ = @SVector [0.01, 0.02, -0.01]
u̇₂ = @SVector [0.013, -0.02, -0.015]

ẇ₁ = @SVector [-0.02423, 0.009047, 0.005896]
ẇ₂ = @SVector [-0.002394, 0.008994, 0.005902]

ü₁ = @SVector [0.008, 0.005, 0.001]
ü₂ = @SVector [0.009, 0.0052, 0.0014]

ẅ₁ = @SVector [-0.002423, 0.009047, 0.005896]
ẅ₂ = @SVector [-0.002394, 0.008994, 0.005902]


Rₑ⁰ = local_R⁰(X₁, X₂)
R₁ = rotation_matrix(Θ₁)
R₂ = rotation_matrix(Θ₂)
ΔR₁ = 1. *ID3
ΔR₂ = 1. *ID3
mat = (E = 1., G = 0.1, Jᵨ = Diagonal(@SVector [20, 10, 10]), Aᵨ = 0.01)
geom = (A = 0.01, J = 0.01, I₃₃ = 0.01, I₂₂ = 0.01)
comp = (nᴳ = 3, zᴳ = @SVector[-sqrt(3/5), 0, sqrt(3/5)], ωᴳ = @SVector[5/9, 8/9, 5/9], kₙ = 10., μ = 0.3, εᵗ = 0.1, ηₙ = 1., damping=0.1)

init = (;X₁, X₂, l₀, Rₑ⁰)

simvars = (;mat, geom, comp, init, sdf=nothing)

sdf = Sphere_SDF{Float64}(1., 1., 0, 0, 0)
simvars_contact = (;mat, geom, comp, init, sdf)

function compare(test, correct, disp=true)
    m = maximum(abs.(correct))
    if disp
        println("-------------------------------------------------------------------------")
        display(test)
        println()
        display(correct)
        println()
        display(test - correct)
        println()
    end
    return maximum(abs.(test - correct))/m
end

# Change of variable 
HH = [I(3)*1. zeros(3,3) zeros(3,3) zeros(3,3); zeros(3,3) Tₛ⁻¹(Θ₁) zeros(3,3) zeros(3,3);zeros(3,3) zeros(3,3) I(3)*1. zeros(3,3);zeros(3,3) zeros(3,3) zeros(3,3) Tₛ⁻¹(Θ₂)]

x₀ = [u₁..., Θ₁..., u₂..., Θ₂...]



function finitediff_statics(exact=true)

    x₀ = [u₁..., Θ₁..., u₂..., Θ₂...]

    _, _, _, ftest, _, _, Ktest, _, _, _, _ = compute(u₁, u₂, R₁, R₂, ΔR₁, ΔR₂, u̇₁, u̇₂, ẇ₁, ẇ₂, ü₁, ü₂, ẅ₁, ẅ₂, simvars, exact, false)

    fun(x) = compute(x[1:3], x[7:9], rotation_matrix(x[4:6]), rotation_matrix(x[10:12]), ΔR₁, ΔR₂, u̇₁, u̇₂, ẇ₁, ẇ₂, ü₁, ü₂, ẅ₁, ẅ₂, simvars, true, false)

    F = FiniteDiff.finite_difference_gradient(x->fun(HH*x)[1], x₀)
    H = FiniteDiff.finite_difference_jacobian(x->fun(HH*x)[4], x₀)

    return [compare(ftest, F), compare(Ktest, H)]

end

function finitediff_dynamic(exact=true)

    x₀ = [u̇₁..., ẇ₁..., u̇₂..., ẇ₂...]

    _, _, _, _, ftest, _, _, _, Mtest, Ctest, _ = compute(u₁, u₂, R₁, R₂, ΔR₁, ΔR₂, u̇₁, u̇₂, ẇ₁, ẇ₂, ü₁, ü₂, ẅ₁, ẅ₂, simvars, exact, true)
    
    funC(x) = compute(u₁, u₂, R₁, R₂, ΔR₁, ΔR₂, x[1:3], x[7:9], x[4:6], x[10:12], ü₁, ü₂, ẅ₁, ẅ₂, simvars, true, true)

    HC = FiniteDiff.finite_difference_jacobian(x->funC(x)[5], x₀)


    x₀ = [ü₁..., ẅ₁..., ü₂..., ẅ₂...]

    
    funM(x) = compute(u₁, u₂, R₁, R₂, ΔR₁, ΔR₂, u̇₁, u̇₂, ẇ₁, ẇ₂, x[1:3], x[7:9], x[4:6], x[10:12], simvars, true, true)

    HM = FiniteDiff.finite_difference_jacobian(x->funM(x)[5], x₀)

    return [compare(Ctest, HC), compare(Mtest, HM)]

end


function finitediff_contact(exact=true)

    x₀ = [u₁..., Θ₁..., u₂..., Θ₂...]

    # frictionless
    newcomp = (comp..., μ = 0., ηₙ = 0., εᵗ = 0.)
    simvars_contact_frictionless = (;mat, geom, comp=newcomp, init, sdf)
    _, _, _, _, _, ftest, _, _, _, _, _ = compute(u₁, u₂, R₁, R₂, ΔR₁, ΔR₂, u̇₁, u̇₂, ẇ₁, ẇ₂, ü₁, ü₂, ẅ₁, ẅ₂, simvars_contact_frictionless, exact, true)
    funK_frictionless(x) = compute(x[1:3], x[7:9], rotation_matrix(x[4:6]), rotation_matrix(x[10:12]), ΔR₁, ΔR₂, u̇₁, u̇₂, ẇ₁, ẇ₂, ü₁, ü₂, ẅ₁, ẅ₂, simvars_contact_frictionless, true, true)
    F = FiniteDiff.finite_difference_gradient(x->funK_frictionless(HH*x)[3], x₀)

    # with friction
    _, _, _, _, _, _, _, Ktest, _, _, Ctest = compute(u₁, u₂, R₁, R₂, ΔR₁, ΔR₂, u̇₁, u̇₂, ẇ₁, ẇ₂, ü₁, ü₂, ẅ₁, ẅ₂, simvars_contact, exact, true)
    funK(x) = compute(x[1:3], x[7:9], rotation_matrix(x[4:6]), rotation_matrix(x[10:12]), ΔR₁, ΔR₂, u̇₁, u̇₂, ẇ₁, ẇ₂, ü₁, ü₂, ẅ₁, ẅ₂, simvars_contact, true, true)
    H = FiniteDiff.finite_difference_jacobian(x->funK(HH*x)[6], x₀)


    x₀ = [u̇₁..., ẇ₁..., u̇₂..., ẇ₂...]

    funC(x) = compute(u₁, u₂, R₁, R₂, ΔR₁, ΔR₂, x[1:3], x[7:9], x[4:6], x[10:12], ü₁, ü₂, ẅ₁, ẅ₂, simvars_contact, true, true)

    HC = FiniteDiff.finite_difference_jacobian(x->funC(x)[6], x₀)

    return [compare(ftest, F), compare(Ktest, H), compare(Ctest, HC)]

end


@test finitediff_statics(true) ≈ [0, 0] atol=1e-7

@test finitediff_dynamic(true) ≈ [0, 0] atol=1e-5

@test finitediff_contact(true) ≈ [0, 0, 0] atol=1e-2



@inferred compute(u₁, u₂, R₁, R₂, ΔR₁, ΔR₂, u̇₁, u̇₂, ẇ₁, ẇ₂, ü₁, ü₂, ẅ₁, ẅ₂, simvars_contact, true, true)

@btime compute($u₁, $u₂, $R₁, $R₂, $ΔR₁, $ΔR₂, $u̇₁, $u̇₂, $ẇ₁, $ẇ₂, $ü₁, $ü₂, $ẅ₁, $ẅ₂, $simvars_contact, $true, $true)