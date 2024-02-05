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
using Statistics
using Plots
import ADerrors
import ADerrors: uwreal, uwerr
using Utils

# Function definitions

## C2/C0 correlator definition. p[1] is the correction to the mass

import Utils: nparameters
mutable struct C2Correlator <: UtilsFunc 
    I
end
(s::C2Correlator)(x,p) = -1/2 * (1/(2*s.I))^2 * p[1]^2 * x^2 + p[2] * x + p[3]
nparameters(::C2Correlator) = 3


## Output functions

qwrite(io, obs::Real) = write(io, "$(obs)")
qwrite(io, obs::uwreal)  = write(io, "$(ADerrors.value(obs)),$(ADerrors.err(obs))")
function qwrite(io, obs::FormalSeries.Series{T,N}) where {T,N} 
    for i in 1:N-1
        qwrite(io, obs[i])
        write(io, ",")
    end
    qwrite(io, obs[N])
end
function qwrite(file, vobs::AbstractArray)
    open(file, "w") do io
        for obs in vobs
            qwrite(io, obs)
            write(io, ",")
        end
    end
end
function qwrite(file, message::String)
    open(file, "w") do io
        write(io, message)
    end
end


# Main body

length(ARGS) == 1 || error("Only one argument is expected! (Path to input file)")
isfile(ARGS[1]) || error("Path provided is not a file")

fname = ARGS[1]
dname = dirname(fname)
fdname = basename(dname)

ens = QuantumRotorExperiments.read_ensemble(fname, QuantumRotor)

iT = ens[1].params.iT
I = ens[1].params.I
nconf = length(ens)

configs = zeros(eltype(ens[1].phi), nconf, iT)

for i in 1:nconf
    for t in 1:iT
        configs[i,t] = sin(ens[i].phi[t])
    end
    print("$i \r")
end

Ct = similar(configs, nconf, 30)
meansphi = zeros(eltype(configs),iT)

for i in 1:iT
    meansphi[i] = mean(configs[:,i])
end

Ct .= 0.0
for t in 1:30
    println(t)
    for it in 10:59
        tup = mod1(it+t-1, iT)
        Ct[:,t] .+= (configs[:,it] .- meansphi[it]) .* (configs[:,tup] .- meansphi[tup]) / 50
    end
end

tag = "test"

y = Vector{FormalSeries.Series{uwreal,3}}(undef, 30)
for i in 1:30
    y[i] = uwreal(Ct[:,i], tag)
end

uwerr.(y)

y0 = [y[i].c[1] for i in 1:length(y)]
y1 = [y[i].c[2] for i in 1:length(y)]
y2 = [y[i].c[3] for i in 1:length(y)]


f = C2Correlator(I)

ys = y2 ./ y0
uwerr.(ys)

xs = collect(0:29)

fitp, cse, cs = fit_data(f, xs, ys, [1.0, 0.0, 0.0])
uwerr.(fitp)

pl = plot_fit(xs, ys, f, fitp)
plot!(pl, xs, ADerrors.value.(ys), yerr=ADerrors.err.(ys), seriestype=:scatter)
plot!(pl, x -> -1/2 * (1/(2*f.I))^2 * 1/pi^2 * x^2, 0, 30, label="Analytic", lw=2)
plot!(pl, title="cse=$(round(cse, sigdigits=2)),cs=$(round(cs, sigdigits=2))")

savefig(pl, joinpath(dname,fdname)*"-fit.pdf")


# Mass order 1

E1 = fitp[1]

# Mass order 0

E0 = log(y0[1]/y0[2])
uwerr(E0)

# Action

uws = uwreal(action.(ens), tag)
uwerr(uws)
uws

# Topological susceptibility

uwchi = uwreal(diff_local_susceptibility.(ens), tag)
uwerr(uwchi)
uwchi


# Dump to file

anafile = joinpath(dname,fdname)*"-ana.csv"
logfile = joinpath(dname,fdname)*"-log.log"

qwrite(anafile, [iT, I, uwchi, E0, E1, uws])

lfilecontent = 
""" The analysis in $anafile contains, in order,

T = $(iT)
I = $(I)
uwchi = $(uwchi)
E0 = $(E0)
E1 = $(E1)
uws = $(uws)


Some extra information that does not appear in the file:

tauint S = $(ADerrors.taui(uws[1], tag))
dtauint S = $(ADerrors.dtaui(uws[1], tag))
"""

qwrite(logfile, lfilecontent)
