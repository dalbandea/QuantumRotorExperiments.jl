
import Base:+

function ADerrors.uwreal(obs::Vector{FormalSeries.Series{T, N}}, ID::String) where {T, N}
    uwobs = FormalSeries.Series{ADerrors.uwreal, N}(ntuple(i -> ADerrors.uwreal([real(obs[j].c[i]) for j in 1:length(obs)], ID), N))
    return uwobs
end

function ADerrors.uwerr(obs::FormalSeries.Series{ADerrors.uwreal, N}) where N
    for i in 1:N
        ADerrors.uwerr(obs[i])
    end
end

+(s::uwreal) = s
