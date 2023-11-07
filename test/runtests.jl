using Test
using LFTSampling
using LFTQuantumRotor
using QuantumRotorExperiments
using FormalSeries


@testset verbose = true "FormalSeries tests" begin

    @testset verbose = true "FS HMC" begin
        include("fshmctests.jl")
    end

    @testset verbose = true "FS I/O" begin
        include("fsiotests.jl")
    end

end

