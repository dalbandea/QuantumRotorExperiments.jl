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

end # module QuantumRotorExperiments
