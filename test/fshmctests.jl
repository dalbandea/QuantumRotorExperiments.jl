
I = 5.0
iT = 100

model = QuantumRotor(
                     Float64, 
                     Series{Float64, 4},
                     I = I, 
                     iT = iT, 
                     BC = PeriodicBC, 
                     disc = StAngleDifferenceDiscretization,
                     theta = Series((0.0, 1.0, 0.0, 0.0))
                    )

randomize!(model)

smplr = HMC(integrator = Leapfrog(1.0, 10))
samplerws = LFTSampling.sampler(model, smplr)
LFTSampling.generate_momenta!(model, samplerws)

@testset verbose = true "$(model.params.BC) HMC reversibility" begin
    model_bckp = deepcopy(model)
    LFTSampling.reversibility!(model, samplerws)
    dphi = model.phi .- model_bckp.phi
    for i in 1:length(model.params.theta.c)
        @test isapprox(zero(model.PRC), mapreduce(x -> abs2(x), +, [dphi[j][i] for j in 1:length(dphi)]), atol = 10.0^(-15+2*i))
    end
end

@testset verbose = true "$(model.params.BC) HMC force" begin
    dF = LFTSampling.force_test(model, samplerws, 1e-6)
    for i in 1:length(model.params.theta.c)
        @test isapprox(zero(model.PRC), dF[i], atol = 1e-5)
    end
end

