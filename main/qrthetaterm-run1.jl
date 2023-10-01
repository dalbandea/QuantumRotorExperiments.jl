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
                     Series{ComplexF64, 4},
                     I = I, 
                     iT = iT, 
                     BC = OpenBC, 
                     disc = StAngleDifferenceDiscretization,
                     theta = Series((0.0im, 1.0+0.0im, 0.0im, 0.0im))
                    )

randomize!(model)

smplr = HMC(integrator = Leapfrog(1.0, 1000))
samplerws = LFTSampling.sampler(model, smplr)

fname = "cfgs-run-I$I-T$iT.bdio"
fb = BDIO_open(fname, "w","Quantum rotor")
BDIO_close!(fb)

@time sample!(model, samplerws)

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



fb2 = BDIO_open(fname, "r")

@time while BDIO_seek!(fb2)
    if BDIO_get_uinfo(fb2) == 1
        BDIO_read(fb2, model2.phi)
    end
end


fb2 = BDIO_open("test-n.bdio", "r")
@time while BDIO_seek!(fb2)
    if BDIO_get_uinfo(fb2) == 1
        res += 1
    end
end
BDIO_close!(fb2)


fb2 = BDIO_open("test.bdio", "r")

BDIO_seek!(fb2)

BDIO_seek!(fb2)

BDIO_read(fb2, model2.phi)

BDIO_close!(fb2)


