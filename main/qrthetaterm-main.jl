# Quantum Rotor
import Pkg
Pkg.activate(".")
using LFTQuantumRotor
using LFTSampling
using QuantumRotorExperiments
using FormalSeries
using BDIO
using ArgParse


# -T 500 -I 5.0 --BC OpenBC --nsteps 500 -n 100000 --wdir /home/david/scratch/projects/phd/8-Strong-CP-Quantum-Rotor/
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "-T"
        help = "lattice size"
        required = true
        arg_type = Int

        "-I"
        help = "I"
        required = true
        arg_type = Float64

        "--BC"
        help = "Boundary conditions: PeriodicBC or OpenBC"
        required = true
        arg_type = String

        "-N"
        help = "theta Taylor order"
        required = true
        arg_type = Int

        "-t"
        help = "tau"
        required = false
        arg_type = Float64
        default = 1.0

        "--nsteps"
        help = "integration steps per trajectory"
        required = true
        arg_type = Int

        "-n"
        help = "number of trajectories"
        required = true
        arg_type = Int

        "--wdir"
        help = "path to directory to save configurations and logs"
        required = true
        arg_type = String
        # default = "configs/"
    end

    return parse_args(s)
end

########################
# PARSE ARGUMENTS ######
########################

parsed_args = parse_commandline()

# Theory Parameters
iT   = parsed_args["T"]
I    = parsed_args["I"]
N    = parsed_args["N"]
BC      = eval(Meta.parse(parsed_args["BC"]))

# HMC parameters
tau     = parsed_args["t"]
nsteps  = parsed_args["nsteps"]
epsilon = tau/nsteps
n_traj  = parsed_args["n"]

# Working directory
wdir    = parsed_args["wdir"]


############
# HMC ######
############

model = QuantumRotor(
                     Float64, 
                     Series{Float64, 3},
                     I = I, 
                     iT = iT, 
                     BC = PeriodicBC, 
                     disc = StAngleDifferenceDiscretization,
                     theta = Series((0.0, 1.0, 0.0))
                    )

randomize!(model)

# Create sampler workspace
smplr = HMC(integrator = Leapfrog(tau, nsteps))
samplerws = LFTSampling.sampler(model, smplr)

fname = joinpath(wdir, "run-QR-FS-theta-N$N-I$I-T$iT-BC$BC-nc$n_traj.bdio")

if isfile(fname)
    error("File already exists. Aborting.")
end

# Thermalize
for i in 1:10000
    @time sample!(model, samplerws)
end

# Run
for i in 1:n_traj
    @time sample!(model, samplerws)
    LFTSampling.save_cnfg(fname, model)
end
