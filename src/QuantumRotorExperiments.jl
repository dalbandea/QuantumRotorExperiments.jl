module QuantumRotorExperiments

import LFTSampling
import LFTSampling: metropolis_accept_reject!
import LFTQuantumRotor
import BDIO
import BDIO: BDIO_write!, BDIO_read
import FormalSeries

include("qrio.jl")
# export BDIO_write!, BDIO_read

function LFTSampling.metropolis_accept_reject!(lftws::L, lftcp::L, dS::FormalSeries.Series) where {L <: LFTQuantumRotor.QuantumRotor} 
    println(dS)
end

end # module QuantumRotorExperiments
