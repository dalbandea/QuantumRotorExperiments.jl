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

function read_ensemble(fname::String, LFT::Type{QuantumRotor}, n::Int64 = 0)
    return LFTSampling.read_ensemble(fname, LFT, n; modul = QuantumRotorExperiments)
end
