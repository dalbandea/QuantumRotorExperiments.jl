
compute_q2(qrws::LFTQuantumRotor.QuantumRotor) = compute_q2(qrws, qrws.params.disc)
function compute_q2(qrws::LFTQuantumRotor.QuantumRotor, disc::Type{D}) where D <: LFTQuantumRotor.AbstractAngleDifferenceDiscretization
    return @views mapreduce(x -> (-im*log(exp(im*x)))^2, +, qrws.phi[1:end-1])/(qrws.params.iT-1)/4pi^2
end

compute_q(qrws::LFTQuantumRotor.QuantumRotor) = compute_q(qrws, qrws.params.disc)
function compute_q(qrws::LFTQuantumRotor.QuantumRotor, disc::Type{D}) where D <: LFTQuantumRotor.AbstractAngleDifferenceDiscretization
    return @views mapreduce(x -> (-im*log(exp(im*x))), +, qrws.phi[1:end-1])/(qrws.params.iT-1)/4pi^2
end

diff_top_charge(qrws::LFTQuantumRotor.QuantumRotor) = diff_top_charge(qrws, qrws.params.disc, qrws.params.BC)
function diff_top_charge(qrws::LFTQuantumRotor.QuantumRotor, disc::Type{D}, ::Type{BC}) where {D <: LFTQuantumRotor.AbstractAngleDifferenceDiscretization, BC <: LFTQuantumRotor.AbstractBoundaryCondition}
    Q = zero(eltype(qrws.phi))
    qt = zero(eltype(qrws.phi))
    for i in 1:qrws.params.iT-1
        qt -= qrws.phi[i]
        Q += sin(qrws.phi[i])
    end

    if BC == LFTQuantumRotor.PeriodicBC
        Q += sin(qt) 
    elseif BC == LFTQuantumRotor.OpenBC
    else
        error("This should not happen")
    end
    return Q/2pi
end
