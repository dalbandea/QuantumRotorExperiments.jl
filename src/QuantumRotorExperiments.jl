module QuantumRotorExperiments

import LFTSampling
import LFTSampling: metropolis_accept_reject!
import LFTQuantumRotor
import BDIO
import BDIO: BDIO_write!, BDIO_read
using FormalSeries
import ADerrors
import ADerrors: uwreal
import LinearAlgebra

include("qrio.jl")
# export BDIO_write!, BDIO_read
include("qrobservables.jl")
export compute_q2, compute_q, diff_top_charge, ZeroMode
include("qraderrors.jl")
include("trivmap.jl")

include("qrseries.jl")
export restart_derivatives, restart_derivatives!
export Mod, Mod!
export diff_local_susceptibility


function LFTSampling.metropolis_accept_reject!(lftws::L, lftcp::L, dS::FormalSeries.Series) where {L <: LFTQuantumRotor.QuantumRotor} 
    # println(dS)
end

generate_pseudofermions!(qrws::QuantumRotor, hmcws::QuantumRotorHMC) = generate_pseudofermions!(qrws,hmcws,qrws.params.disc,qrws.params.BC,qrws.aux)
function generate_pseudofermions!(qrws::QuantumRotor, hmcws::QuantumRotorHMC, disc::Type{D}, BC::Type{B}, aux::AUX) where {D <: LFTQuantumRotor.AbstractDiscretization, B <: LFTQuantumRotor.AbstractBoundaryCondition, AUX <: AbstractAuxFields} end


end # module QuantumRotorExperiments
