using DifferentialEquations
include("swave.jl")

#Define simulation parameters
N = 50000; #number of spins
D = 10.0 #Half the energy bandwidth
δ = 2*D/(N-1) #Level spacing
ν = 1/δ

t0 = 0.0; #start time
tf = 100.0; #final time
tspan = (t0,tf);

Δ0 = .6 #1.0*exp(-π/2)*(1-.1)
Δf = .8 #1.0

ε = -D .+ δ*Vector(0:N-1);

g0 = swave.get_g(Δ0,N,ε)
gf = swave.get_g(Δf,N,ε)

#Create rank 2 tensor holding system state
S = Array{Float64,2}(undef,N,3)

for i in 1:N
    denom = 2*sqrt(ε[i]^2+Δ0^2)
    S[i,1] = Δ0/denom
    S[i,2] = 0.0
    S[i,3] = -ε[i]/denom
end

# Print parameters
println("#N: ",N)
println("#D0: ",Δ0)
println("#Df: ",Δf)
println("#g0: ",g0)
println("#gf: ",gf)
print("#")

p = swave.params(N,g0,Δ0,gf,Δf,ε);
prob = ODEProblem(swave.eom,S,tspan,p);
Δ_t = SavedValues(Float64, Tuple{Float64,Float64});
cb = SavingCallback((u,t,integrator)->(swave.calcΔ(u,gf,p)), Δ_t, saveat=0.0:0.01:tf);
@time sol = solve(prob,DP5(),save_everystep=false,callback=cb,maxiters=1e10,reltol=1e-4,abstol=1e-9);

Dx = Array{Float64,1}(undef,length(Δ_t.t))
Dy = Array{Float64,1}(undef,length(Δ_t.t))
for i in 1:length(Δ_t.t)
    println(Δ_t.t[i], ",", Δ_t.saveval[i][1])
end
