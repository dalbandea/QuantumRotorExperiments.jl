import FormalSeries: Series
import Base: mod
import LFTQuantumRotor: QuantumRotor, Mod
import Base: zero
import Base: broadcastable
import FormalSeries: genseries

function restart_derivatives(s::FormalSeries.Series{T,N}) where {T, N}
    return FormalSeries.Series{T, N}(ntuple(i -> i == 1 ? s[1] : zero(T), N))
end

function restart_derivatives!(qrws::QuantumRotor)
    qrws.phi .= restart_derivatives.(qrws.phi)
    return nothing
end

function Mod(s::FormalSeries.Series{T,N}, alpha) where {T, N}
    return FormalSeries.Series{T, N}(ntuple(i -> i == 1 ? Mod(s[1],alpha) : s[i], N))
end

function Mod!(qrws::QuantumRotor, alpha)
    qrws.phi .= Mod.(qrws.phi, alpha)
    return nothing
end

function mod(z::T, alpha) where T <: Complex
    imag(z) == zero(imag(z)) || error("Non-zero imaginary part in mod")
    return mod(real(z), alpha)
end


Base.:+(s1::FormalSeries.Series{uwreal,N}, s2::uwreal)      where {N} = genseries(FormalSeries.Series{uwreal,N}, i -> i == 1 ? s1.c[1] + s2 : s1.c[i])

Base.zero(::Type{uwreal}) = uwreal([0.0, 0.0], "noerror")

broadcastable(x::FormalSeries.Series) = Ref(x)
