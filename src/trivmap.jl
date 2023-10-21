import LFTQuantumRotor: AbstractAuxFields, FallbackAuxField, StAngleDifferenceDiscretization

struct uQRTrivMap{A1 <: AbstractArray, A2 <: AbstractArray} <: AbstractAuxFields
    C::A1
    psi::A2
end
export uQRTrivMap

##########
# Action #
##########

import LFTQuantumRotor: QuantumRotor, PeriodicBC, QuantumRotorHMC
import LFTQuantumRotor: action, action_t, boundary_action, action_discretization_factor

function action(qrws::QuantumRotor, aux::uQRTrivMap)
    qrws.phi .= aux.C * aux.psi
    return action(qrws, FallbackAuxField())
end


#######
# HMC #
#######

import LFTQuantumRotor: generate_momenta!, force!, boundary_force!
import LFTQuantumRotor: update_fields!
import LFTSampling: generate_pseudofermions!


generate_pseudofermions!(qrws::QuantumRotor, hmcws::QuantumRotorHMC) = generate_pseudofermions!(qrws,hmcws,qrws.params.disc,qrws.params.BC,qrws.aux)
function generate_pseudofermions!(qrws::QuantumRotor, hmcws::QuantumRotorHMC, disc::Type{D}, BC::Type{B}, aux::AUX) where {D <: LFTQuantumRotor.AbstractAngleDifferenceDiscretization, B <: LFTQuantumRotor.AbstractBoundaryCondition, AUX <: AbstractAuxFields} end
function generate_pseudofermions!(qrws::QuantumRotor, hmcws::QuantumRotorHMC, disc::Type{StAngleDifferenceDiscretization}, BC::Type{PeriodicBC}, aux::uQRTrivMap)
    aux.psi .= transpose(aux.C) * qrws.phi
    return nothing
end

function force!(qrws::QuantumRotor, hmcws::QuantumRotorHMC, disc::Type{StAngleDifferenceDiscretization}, BC::Type{PeriodicBC}, aux::uQRTrivMap)

    qrws.phi .= aux.C * aux.psi
    sumphi = @views sum(qrws.phi[1:end-1])
    hmcws.frc .= -qrws.params.I * (transpose(aux.C) * sin.(qrws.phi) .+ 
                                   sin(sumphi) * transpose(reduce(+, aux.C, dims=1)))
    hmcws.frc[end] = zero(eltype(qrws.phi))
    return nothing
end

function update_fields!(qrws::QuantumRotor, epsilon, hmcws::QuantumRotorHMC, aux::uQRTrivMap)
    # Update psi field
    # aux.psi .= transpose(aux.C) * qrws.phi
    aux.psi .= aux.psi .+ epsilon .* hmcws.mom ./ hmcws.params.width^2
    return nothing
end


####################
# Trivializing map #
####################


function orthonormalizing_base(T)
    S = zeros(T-1,T-1)
    res = zeros(T,T)
    S[1,:] .= -1.0
    S[:,1] .= 1.0
    for i in 2:T-1
        S[i,i] = 1.0
    end
    res[1:end-1,1:end-1] .= copy(LinearAlgebra.qr(S).Q)
    return res
end
export orthonormalizing_base
