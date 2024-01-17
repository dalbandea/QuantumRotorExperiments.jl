# Quantum Rotor
using Revise
import Pkg
Pkg.activate(".")
using LFTQuantumRotor
using LFTSampling
using QuantumRotorExperiments
import FormalSeries
using BDIO
using ProgressBars

I = 5.0
iT = 100

model = QuantumRotor(
                     Float64, 
                     FormalSeries.Series{Float64, 3},
                     I = I, 
                     iT = iT, 
                     BC = OpenBC, 
                     disc = StandardDiscretization,
                     theta = FormalSeries.Series((0.0,1.0, 0.0)),
                    )

randomize!(model)

smplr = HMC(integrator = OMF4(1.0, 15))
samplerws = LFTSampling.sampler(model, smplr)

for i in 1:1000
    @time sample!(model, samplerws)
end

restart_derivatives!(model)
Mod!(model, 2pi)

ens = [deepcopy(model) for i in 1:1000000]

for i in ProgressBar(1:1000000)
    for j in 1:10
        sample!(model, samplerws)
    end
    sample!(model, samplerws)
    Mod!(model, 2pi)
    ens[i].phi .= model.phi
    # print("$i / 100000 \r")
end

save_ensemble("BIG-SIM-I5-T100.bdio", ens)
