module swave

struct params
    N::Int64
    g0::Float64
    Δ0::ComplexF64
    gf::Float64
    Δf::ComplexF64
    ε::Array{Float64,1}
end

function calcΔ(S::Array{Float64,2},g::Float64,p::params)
    Δxsum = 0; Δysum = 0;
    @inbounds for i in 1:p.N
        Δxsum += S[i,1]
        Δysum += S[i,2]
    end
    return (g*Δxsum,g*Δysum)
end

function eom(du::Array{Float64,2},u::Array{Float64,2},p::params,t::Float64)

    Δx,Δy = calcΔ(u,p.gf,p)

    @inbounds for i in 1:p.N
        du[i,1] = -2*(p.ε[i]*u[i,2] + u[i,3]*Δy)
        du[i,2] =  2*(p.ε[i]*u[i,1] + u[i,3]*Δx)
        du[i,3] = -2*(Δx*u[i,2] - Δy*u[i,1])
    end

end

function get_g(Δ0::Float64,N::Int64,ε::Array{Float64,1})
    res = 0.0
    @inbounds for i in 1:N
        res += 1.0 / sqrt(ε[i]^2 + Δ0^2)
    end

    g = 2.0 / res

    return g
end
end
