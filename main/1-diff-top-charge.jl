# Quantum Rotor
using Revise
import Pkg
Pkg.activate(".")
using LFTQuantumRotor
using LFTSampling
using QuantumRotorExperiments
using FormalSeries
using BDIO

I = 5.0
iT = 500

model = QuantumRotor(
                     Float64, 
                     # Series{ComplexF64, 4},
                     I = I, 
                     iT = iT, 
                     BC = PeriodicBC,
                     disc = StAngleDifferenceDiscretization,
                     # theta = Series((0.0im, 1.0+0.0im, 0.0im, 0.0im))
                    )

randomize!(model)

smplr = HMC(integrator = Leapfrog(1.0, 30))
samplerws = LFTSampling.sampler(model, smplr)

fname = "cfgs-run-I$I-T$iT.bdio"
fb = BDIO_open(fname, "w","Quantum rotor")
BDIO_close!(fb)

@time sample!(model, samplerws, do_winding=true)

Ss = Vector{Float64}(undef, 100000)
Qs = Vector{Float64}(undef, 100000)
Qs2 = Vector{Float64}(undef, 100000)

for i in 1:100000
    @time sample!(model, samplerws, do_winding=true)
    Ss[i] = action(model)
    Qs[i] = LFTQuantumRotor.top_charge(model)
    Qs2[i] = diff_top_charge(model)
end

using ADerrors

id = "test"

uws = uwreal(Ss, id)
uwerr(uws)
uws
taui(uws,id)
dtaui(uws,id)

cosphi = 1-uws/(I*(iT))
uwerr(cosphi)
cosphi

uwchi = uwreal(Qs2.^2/iT, id)
uwerr(uwchi)
uwchi

tau = taui(cosphi, id)
dtau = dtaui(cosphi, id)



for i in 1:1000000
    @time begin
        sample!(model, samplerws)
        fb = BDIO_open(fname, "a","Quantum rotor")
        BDIO_start_record!(fb, BDIO_BIN_F64LE, 1, true)
        BDIO_write!(fb,model.phi)
        BDIO_write_hash!(fb)
        BDIO_close!(fb)
    end
end


for i in 1:1000000
    @time begin
        sample!(model, samplerws)
        fb = BDIO_open(fname, "a","Quantum rotor")
        BDIO_start_record!(fb, BDIO_BIN_F64LE, 1, true)
        BDIO_write!(fb,model.phi)
        BDIO_write_hash!(fb)
        BDIO_close!(fb)
    end
end

