# Quantum Rotor
using Revise
import Pkg
Pkg.activate(".")
using LFTQuantumRotor, LFTSampling, QuantumRotorExperiments
using FormalSeries
using BDIO
using ADerrors

I = 5.0
iT = 500

model = QuantumRotor(
                     Float64, 
                     Series{ComplexF64, 4},
                     I = I, 
                     iT = iT, 
                     BC = OpenBC, 
                     disc = StAngleDifferenceDiscretization,
                     theta = Series((0.0im, 1.0+0.0im, 0.0im, 0.0im))
                    )

model.phi .= zero(model.PRC)

fname = "cfgs-run-I$I-T$iT.bdio"

fb2 = BDIO_open(fname, "r")

q2s = Vector{Series{ComplexF64,4}}()
qs = Vector{Series{ComplexF64,4}}()

while BDIO_seek!(fb2)
    if BDIO_get_uinfo(fb2) == 1
        BDIO_read(fb2, model.phi)
        push!(q2s, compute_q2(model))
        push!(qs, compute_q(model))
    end
end

BDIO_close!(fb2)

ID = "test"

uwq2s = uwreal(q2s, ID)

uwerr(uwq2s)

uwqs = uwreal(qs, ID)

uwerr(uwqs)

uwchi = uwq2s - uwqs*uwqs

uwerr(uwchi)

