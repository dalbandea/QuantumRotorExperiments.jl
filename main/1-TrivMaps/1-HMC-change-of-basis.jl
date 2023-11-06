# Quantum Rotor
using Revise
import Pkg
Pkg.activate(".")
using LFTQuantumRotor
using LFTSampling
using QuantumRotorExperiments
using FormalSeries
using BDIO

I = 10.0
iT = 1000

C = orthonormalizing_base(iT)
psi = zeros(iT)
auxtool = uQRTrivMap(C, psi)

model = QuantumRotor(
                     Float64, 
                     aux = auxtool,
                     I = I, 
                     iT = iT, 
                     BC = PeriodicBC,
                     disc = StAngleDifferenceDiscretization,
                    )

randomize!(model)

smplr = HMC(integrator = Leapfrog(1.0, 70))
samplerws = LFTSampling.sampler(model, smplr)


Ss = Vector{Float64}(undef, 1000000)
Qs = Vector{Float64}(undef, 1000000)

@time sample!(model, samplerws, do_winding=false)


LFTQuantumRotor.top_charge(model)

for i in 1:1000000
    @time sample!(model, samplerws, do_winding=false)
    # Ss[i] = action(model)
    Ss[i] = action(model, LFTQuantumRotor.FallbackAuxField())
    Qs[i] = LFTQuantumRotor.top_charge(model)
    println(i)
end


using ADerrors

id = "test2"

uws = uwreal(Ss, id)
uwcos = 1 - uws/(iT*I)
uwerr(uwcos)
uwcos

uwchi = uwreal(Qs.^2/iT, id)
uwerr(uwchi)
uwchi

dtaui(uwchi, id)


# Force test

randomize!(model)

s1 = action(model)
LFTQuantumRotor.force!(model, samplerws)
epsi = 0.0001
model.phi[2] += epsi
s2 = action(model)

(s2-s1)/epsi
samplerws.frc[2]


# Force test psi

randomize!(model)
model.aux.psi .= transpose(model.aux.C) * model.phi
s1 = action(model)
LFTQuantumRotor.force!(model, samplerws)
epsi = 0.0001
model.aux.psi[2] += epsi
# model.phi .= model.aux.C * model.aux.psi
s2 = action(model)
(s2-s1)/epsi
samplerws.frc[2]

randomize!(model)
model.aux.psi .= transpose(model.aux.C) * model.phi

smplr = HMC(integrator = Leapfrog(1.0, 30))
samplerws = LFTSampling.sampler(model, smplr)

@time sample!(model, samplerws)



randomize!(model)
model.aux.psi .= transpose(model.aux.C) * model.phi

modelbk = deepcopy(model)

LFTSampling.leapfrog!(model, samplerws, 0.01, 1000)

samplerws.mom .*= -1

LFTSampling.leapfrog!(model, samplerws, 0.01, 1000)

model.phi .- modelbk.phi

