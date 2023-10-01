# Quantum Rotor
using Revise
import Pkg
Pkg.activate(".")
using LFTQuantumRotor
using LFTSampling
using QuantumRotorExperiments
using FormalSeries
using BDIO

I = 1.0
iT = 100

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

model.phi .= real.(model.phi)

action(model)

LFTQuantumRotor.top_charge(model)

smplr = HMC(integrator = Leapfrog(1.0, 1000))
samplerws = LFTSampling.sampler(model, smplr)


Xs = []
for i in 1:10000
    @time sample!(model, samplerws)
    push!(Xs, mapreduce(x->x^2,+,model.phi[1:end-1])/(model.params.iT-1))
end


for i in 1:model.params.iT
    model.phi[i] = real(model.phi[i][1])
end

@time sample!(model, samplerws)

model.phi .= -im*log.(exp.(im*model.phi))

model.phi

model2 = deepcopy(model)
model2.phi .= 0.0

fb = BDIO_open("test.bdio", "a","Quantum rotor")
BDIO_start_record!(fb, BDIO_BIN_F64LE, 1, true)
BDIO_write!(fb,model.phi)
BDIO_write_hash!(fb)
BDIO_close!(fb)


fb2 = BDIO_open("test.bdio", "r")

@time while BDIO_seek!(fb2)
    if BDIO_get_uinfo(fb2) == 1
        BDIO_read(fb2, model2.phi)
    end
end

BDIO_close!(fb2)


fb2 = BDIO_open("test.bdio", "r")

BDIO_seek!(fb2)

BDIO_seek!(fb2)

BDIO_read(fb2, model2.phi)

BDIO_close!(fb2)

# HMC OBC CP {{{

I = 5.0
iT = 500
model = QuantumRotor(
                     Float64, 
                     I = I, 
                     iT = iT, 
                     BC = OpenBC, 
                     disc = CPAngleDifferenceDiscretization,
                     # disc = StAngleDifferenceDiscretization,
                    )


model.phi .= 0.0

smplr = HMC(integrator = Leapfrog(1.0, 15))
samplerws = LFTSampling.sampler(model, smplr)

sample!(model, samplerws)


Ss = Vector{Float64}()

for i in 1:100000
    @time sample!(model, samplerws)
    push!(Ss, mapreduce(x -> cos(x), +, model.phi[1:end-1])/(iT-1))
end

using ADerrors

id = "test"
uws = uwreal(Ss, id)
uwerr(uws)
uws

cosphi = 1-uws/(I*(iT-1))

# }}}
