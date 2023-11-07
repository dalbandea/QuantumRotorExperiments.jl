import BDIO: BDIO_write!, BDIO_read
import LFTSampling: save_cnfg_header, read_cnfg_info
import LFTQuantumRotor: QuantumRotor
import FormalSeries: Series

function BDIO.BDIO_write!(fb::BDIO.BDIOstream, vec::Vector{FormalSeries.Series{T,N}}) where {T <: Real, N}
    for i in 1:length(vec), n in 1:N
        BDIO.BDIO_write!(fb,[vec[i][n]])
    end
end

function BDIO.BDIO_read(fb::BDIO.BDIOstream, vec::Vector{FormalSeries.Series{T,N}}) where {T <: Real, N}
    reg = zeros(Float64,N)
    for i in 1:length(vec)
        BDIO.BDIO_read(fb,reg)
        vec[i] = FormalSeries.Series{T,N}(ntuple(n -> reg[n], N))
    end
end

function BDIO.BDIO_write!(fb::BDIO.BDIOstream, vec::Vector{FormalSeries.Series{T,N}}) where {T <: Complex, N}
    for i in 1:length(vec), n in 1:N
        BDIO.BDIO_write!(fb,[real(vec[i][n])])
        BDIO.BDIO_write!(fb,[imag(vec[i][n])])
    end
end

function BDIO.BDIO_read(fb::BDIO.BDIOstream, vec::Vector{FormalSeries.Series{T,N}}) where {T <: Complex, N}
    reg = zeros(Float64,2*N)
    for i in 1:length(vec)
        BDIO.BDIO_read(fb,reg)
        vec[i] = FormalSeries.Series{T,N}(ntuple(n -> complex(reg[2*n-1], reg[2*n]), N))
    end
end

BDIO.BDIO_write!(fb::BDIO.BDIOstream, qrws::LFTQuantumRotor.QuantumRotor) = BDIO.BDIO_write!(fb, qrws.phi)

BDIO.BDIO_read(fb::BDIO.BDIOstream, qrws::LFTQuantumRotor.QuantumRotor) = BDIO.BDIO_read(fb, qrws.phi)

function save_cnfg_header(fb::BDIO.BDIOstream, qrws::QuantumRotor)
    BDIO.BDIO_write!(fb, [qrws.params.I])
    BDIO.BDIO_write!(fb, [convert(Int32, qrws.params.iT)])
    BDIO_write!(fb, string(qrws.params.BC)*"\0")
    BDIO_write!(fb, string(typeof(qrws.params.theta))*"\0")
    BDIO.BDIO_write!(fb, [qrws.params.theta])
    return nothing
end

function read_cnfg_info(fname::String, ::Type{LFTQuantumRotor.QuantumRotor})

    fb = BDIO.BDIO_open(fname, "r")

    while BDIO.BDIO_get_uinfo(fb) != 1
        BDIO.BDIO_seek!(fb)
    end

    ifoo    = Vector{Float64}(undef, 1)
    BDIO.BDIO_read(fb, ifoo)
    I    = ifoo[1]
    ifoo    = Vector{Int32}(undef, 1)
    BDIO.BDIO_read(fb, ifoo)
    iT      = convert(Int64, ifoo[1])
    BC      = eval(Meta.parse(BDIO.BDIO_read_str(fb)))
    thtp    = eval(Meta.parse(BDIO.BDIO_read_str(fb)))
    thfoo   = [zero(thtp)]
    BDIO.BDIO_read(fb, thfoo)
    theta   = thfoo[1]
    println(BC)

    model = LFTQuantumRotor.QuantumRotor(
                         Float64, 
                         thtp,
                         I = I, 
                         iT = iT, 
                         BC = BC, 
                         disc = StAngleDifferenceDiscretization,
                         theta = theta
                        )

    return fb, model
end
