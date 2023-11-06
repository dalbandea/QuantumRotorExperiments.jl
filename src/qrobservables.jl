
compute_q2(qrws::LFTQuantumRotor.QuantumRotor) = compute_q2(qrws, qrws.params.disc)
function compute_q2(qrws::LFTQuantumRotor.QuantumRotor, disc::Type{D}) where D <: LFTQuantumRotor.AbstractAngleDifferenceDiscretization
    return @views mapreduce(x -> (-im*log(exp(im*x)))^2, +, qrws.phi[1:end-1])/(qrws.params.iT-1)/4pi^2
end

compute_q(qrws::LFTQuantumRotor.QuantumRotor) = compute_q(qrws, qrws.params.disc)
function compute_q(qrws::LFTQuantumRotor.QuantumRotor, disc::Type{D}) where D <: LFTQuantumRotor.AbstractAngleDifferenceDiscretization
    return @views mapreduce(x -> (-im*log(exp(im*x))), +, qrws.phi[1:end-1])/(qrws.params.iT-1)/4pi^2
end


