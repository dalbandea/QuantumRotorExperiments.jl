
compute_q2(qrws::LFTQuantumRotor.QuantumRotor) = compute_q2(qrws, qrws.params.disc)
function compute_q2(qrws::LFTQuantumRotor.QuantumRotor, disc::Type{D}) where D <: LFTQuantumRotor.AbstractAngleDifferenceDiscretization
    return @views mapreduce(x -> (-im*log(exp(im*x)))^2, +, qrws.phi[1:end-1])/(qrws.params.iT-1)/4pi^2
end

compute_q(qrws::LFTQuantumRotor.QuantumRotor) = compute_q(qrws, qrws.params.disc)
function compute_q(qrws::LFTQuantumRotor.QuantumRotor, disc::Type{D}) where D <: LFTQuantumRotor.AbstractAngleDifferenceDiscretization
    return @views mapreduce(x -> (-im*log(exp(im*x))), +, qrws.phi[1:end-1])/(qrws.params.iT-1)/4pi^2
end



diff_local_susceptibility(qrws::LFTQuantumRotor.QuantumRotor) = diff_local_susceptibility(qrws, qrws.params.disc, qrws.params.BC)
function diff_local_susceptibility(qrws::LFTQuantumRotor.QuantumRotor, disc::Type{LFTQuantumRotor.StandardDiscretization}, ::Type{BC}) where BC <: LFTQuantumRotor.AbstractBoundaryCondition
    chi = zero(eltype(qrws.phi))
    for i in 1:qrws.params.iT-1
        chi += sin(qrws.phi[i+1]-qrws.phi[i])^2
    end

    if BC == LFTQuantumRotor.PeriodicBC
        chi += sin(qrws.phi[1]-qrws.phi[qrws.params.iT])^2
    end

    return chi/(2pi)^2/(qrws.params.iT-1)
end

import LFTQuantumRotor: AbstractAuxFields

struct ZeroMode <: AbstractAuxFields end

import LFTQuantumRotor: left_boundary_force!, update_fields!

function left_boundary_force!(qrws::QuantumRotor, hmcws::LFTQuantumRotor.QuantumRotorHMC, disc::Type{LFTQuantumRotor.StandardDiscretization}, BC::Type{LFTQuantumRotor.OpenBC}, aux::ZeroMode)
    hmcws.frc[1] = 0.0
    println("hi")
    return nothing
end

function update_fields!(qrws::QuantumRotor, epsilon, hmcws::LFTQuantumRotor.QuantumRotorHMC, aux::ZeroMode)
    # Update phi field
    qrws.phi .= qrws.phi .+ epsilon .* hmcws.mom ./ hmcws.params.width^2
    qrws.phi[1] = 0.0
    return nothing
end
