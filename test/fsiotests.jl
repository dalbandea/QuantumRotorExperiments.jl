
I = 5.0
iT = 100

model = QuantumRotor(
                     Float64, 
                     Series{Float64, 3},
                     I = I, 
                     iT = iT, 
                     BC = PeriodicBC, 
                     disc = StAngleDifferenceDiscretization,
                     theta = Series((0.0, 1.0, 0.0))
                    )

LFTQuantumRotor.randomize!(model)
fname = "QRiotest.bdio"
isfile(fname) && error("File already exists!")
LFTSampling.save_cnfg(fname, model)
fb, model2 = QuantumRotorExperiments.read_cnfg_info(fname, QuantumRotor)
LFTQuantumRotor.read_next_cnfg(fb, model2)
rm(fname, force=true)

@testset "Periodic BC I/O" begin
    @test model.params == model2.params
    @test model.phi == model2.phi
end
