# Quantum Rotor
import Pkg
Pkg.activate(".")
using LFTSampling
using LFTQuantumRotor
using QuantumRotorExperiments
using TOML
using Dates
import FormalSeries

infile = "main/infile.in"
pdata = TOML.parsefile(infile)

# Read model parameters

I = pdata["Model params"]["I"]
iT = pdata["Model params"]["T"]
BC = eval(Meta.parse(pdata["Model params"]["BC"]))
disc = eval(Meta.parse(pdata["Model params"]["disc"]))
theta = eval(Meta.parse(pdata["Model params"]["theta"]))

# Read HMC parameters

integrator = eval(Meta.parse(pdata["HMC params"]["integrator"]))
tau = pdata["HMC params"]["tau"]
nsteps = pdata["HMC params"]["nsteps"]
ntherm = pdata["HMC params"]["ntherm"]
ntraj = pdata["HMC params"]["ntraj"]
discard = pdata["HMC params"]["discard"]
windings = pdata["HMC params"]["windings"]

# Working directory

wdir = pdata["Working directory"]["wdir"]

dt = Dates.now()
wdir_sufix = "_D"*Dates.format(dt, "yyyy-mm-dd-HH-MM-SS.ss")
fname = "cfgs-run-I$I-T$iT"*wdir_sufix

fdir = joinpath(wdir, fname)
configfile = joinpath(fdir, fname*".bdio")
mkpath(fdir)
cp(infile, joinpath(fdir,splitpath(infile)[end]))

model = QuantumRotor(
                     Float64, 
                     typeof(theta),
                     I = I, 
                     iT = iT, 
                     BC = BC, 
                     disc = disc,
                     theta = theta
                    )

randomize!(model)

smplr = HMC(integrator = integrator(tau, nsteps))
samplerws = LFTSampling.sampler(model, smplr)

# Thermalization

for i in 1:ntherm
    @time sample!(model, samplerws, do_winding=windings)
    Mod!(model, 2pi)
end


# Run

ens = [deepcopy(model) for i in 1:ntraj]

@time for i in 1:ntraj
    for j in 1:discard
        @time sample!(model, samplerws, do_winding=windings)
    end
    @time sample!(model, samplerws, do_winding=windings)
    Mod!(model, 2pi)
    ens[i].phi .= model.phi
end

save_ensemble(configfile, ens)

