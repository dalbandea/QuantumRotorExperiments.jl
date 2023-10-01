
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
