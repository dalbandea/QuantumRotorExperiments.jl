# Quantum Rotor
using Revise
import Pkg
Pkg.activate(".")
using LFTSampling
using LFTQuantumRotor
using QuantumRotorExperiments
using FormalSeries
using BDIO

## {{{

fname = "/home/david/scratch/projects/phd/8-Strong-CP-Quantum-Rotor/run-QR-FS-theta-N3-I5.0-T500-BCOpenBC-nc500000.bdio"

ens = LFTSampling.read_ensemble(fname, QuantumRotor, 100000)


aconf = similar(ens[1].phi)


function angdiff_to_ang!(aconf, dconf)
    aconf[1] = 0.0
    for i in 2:length(aconf)
        aconf[i] = aconf[i-1] + dconf[i]
    end
end


Ct = similar(aconf, 20000, 30)
Ct .= 0.0

phmean = similar(aconf, 30)
phmean .= 0.0

for i in 1:20000
    angdiff_to_ang!(aconf, ens[i].phi)
    for t in 1:30
        phmean[t] += aconf[t+50] / 20000
    end
    print("$i \r")
end

for i in 1:20000
    angdiff_to_ang!(aconf, ens[i].phi)
    for t in 1:30
        for it in 1:20
            Ct[i, t] += aconf[it+50] * aconf[it+50+t] / 20
        end
    end
    print("$i \r")
end


configs = zeros(Float64, 20000, 500)

for i in 1:20000
    angdiff_to_ang!(aconf, ens[i].phi)
    for t in 1:500
        configs[i,t] = aconf[t][1]
    end
    print("$i \r")
end

configs = LFTQuantumRotor.Mod.(configs, 2pi)


iT = ens[1].params.iT

meansphi = zeros(iT)
for i in 1:iT
    meansphi[i] = mean(configs[:,i])
end


Ct = similar(configs, 20000, 50)
Ct .= 0.0


for t in 1:50
    println(t)
    for it in 1:30
        tup = it+t-1+50
        Ct[:,t] .+= (configs[:,it+50] .- meansphi[it+50]) .* (configs[:,tup+50] .- meansphi[tup+50]) / 30
    end
end

Ctmean = zeros(30)
for i in 1:30
    Ctmean[i] = mean(Ct[:,i])
end


for i in 698:20000
    println(i)
    if sum(isnan.(configs[i,:])) > 0
        break
    end
end

cfg = [ens[803].phi[i][1] for i in 1:500]
for i in 1:500
    println(i)
    if isnan(cfg[i])
        break
    end
end

## }}}


I = 5.0
iT = 100

model = QuantumRotor(
                     Float64, 
                     Series{Float64, 3},
                     I = I, 
                     iT = iT, 
                     BC = OpenBC, 
                     disc = StAngleDifferenceDiscretization,
                     theta = Series((0.0, 1.0, 0.0))
                    )

randomize!(model)

model.phi .= 0.0

smplr = HMC(integrator = Leapfrog(1.0, 20))
samplerws = LFTSampling.sampler(model, smplr)

fname = "cfgs-run-I$I-T$iT.bdio"
fb = BDIO_open(fname, "w","Quantum rotor")
BDIO_close!(fb)

for i in 1:10000
    @time sample!(model, samplerws)
end

ens = [deepcopy(model) for i in 1:10000]

for i in 1:10000
    @time sample!(model, samplerws)
    ens[i].phi .= model.phi
    if any(isnan.(action(model).c))
        println(i)
        break
    end
end

findmax([action(ens[i])[3] for i in 1:10000])



function angdiff_to_ang!(aconf, dconf)
    aconf[1] = 0.0
    for i in 2:length(aconf)
        aconf[i] = aconf[i-1] + dconf[i]
    end
end

aconf = similar(ens[1].phi)

configs = zeros(Float64, 100000, 100)

for i in 1:100000
    angdiff_to_ang!(aconf, ens[i].phi)
    for t in 1:100
        configs[i,t] = aconf[t][1]
    end
    print("$i \r")
end

configs = LFTQuantumRotor.Mod.(configs, 2pi)


iT = ens[1].params.iT

meansphi = zeros(iT)
for i in 1:iT
    meansphi[i] = mean(configs[:,i])
end


Ct = similar(configs, 10000, 40)
Ct .= 0.0


for t in 1:40
    println(t)
    for it in 20:50
        tup = it+t-1
        Ct[:,t] .+= (configs[:,it] .- meansphi[it]) .* (configs[:,tup] .- meansphi[tup]) / 30
    end
end

Ctmean = zeros(40)
for i in 1:40
    Ctmean[i] = mean(Ct[:,i])
end


using ADerrors

y = Vector{uwreal}(undef, )

using Plots

plot(Ctmean)

Ms = similar(Ctmean, 39)
for i in 1:length(Ms)
    Ms[i] = log(Ctmean[i]/Ctmean[i+1])
end

plot(Ms)


# With FormalSeries

ens_bck = deepcopy(ens)
ens_ang = deepcopy(ens)

aconf = similar(ens[1].phi)

configs = zeros(Float64, 10000, 100)

for i in 1:10000
    angdiff_to_ang!(aconf, ens_bck[i].phi)
    ens_ang[i].phi .= aconf
    print("$i \r")
end

for i in 1:length(ens)
    for t in 1:100
        ens[i].phi[t] = Series{Float64, 3}(ntuple(n -> n == 1 ? LFTQuantumRotor.Mod(ens_ang[i].phi[t][1],2pi) : ens_ang[i].phi[t][n], 3))
    end
    println(i)
end


meansphi = similar(ens[1].phi, iT)
meansphi .= 0.0

for i in 1:10000
    for t in 1:iT
        meansphi[t] += ens[i].phi[t] / 10000
    end
    print("$i \r")
end

Ct = similar(ens[1].phi, 10000, 10)
Ct .= 0.0


for i in 1:10000
    for t in 1:10
        for it in 40:50
            tup = it+t-1
            Ct[i,t] += (ens[i].phi[it] - meansphi[it]) * (ens[i].phi[tup] - meansphi[tup]) / 11
        end
    end
    println(i)
end


Ctmean = similar(Ct,10)
for i in 1:40
    Ctmean[i] = mean(Ct[:,i])
end

plot([el[3] for el in Ctmean])
