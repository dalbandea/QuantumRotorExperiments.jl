module QuantumRotorExperiments

import LFTSampling
import LFTSampling: metropolis_accept_reject!
import LFTQuantumRotor
import BDIO
import BDIO: BDIO_write!, BDIO_read
import FormalSeries
import ADerrors
import ADerrors: uwreal

include("qrio.jl")
# export BDIO_write!, BDIO_read
include("qrobservables.jl")
export compute_q2, compute_q, diff_top_charge
include("qraderrors.jl")


function LFTSampling.metropolis_accept_reject!(lftws::L, lftcp::L, dS::FormalSeries.Series) where {L <: LFTQuantumRotor.QuantumRotor} 
    println(dS)
end

end # module QuantumRotorExperiments
