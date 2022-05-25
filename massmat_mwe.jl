using OrdinaryDiffEq
using LinearSolve
using ComponentArrays

function f(resid,du,u,p,t)
    resid[1] = - 0.04*u[1] + 1e4*u[2]*u[3] - du[1]
    resid[2] = + 0.04*u[1] - 3e7*u[2]^2 - 1e4*u[2]*u[3] - du[2]
    resid[3] = u[1] + u[2] + u[3] - 1.0
end

u₀ = [1.0,0.0,0.0]
du₀ = [0.0,0.0,0.0]


prob = DAEProblem(f,du₀,u₀,tspan,differential_vars=[true,true,false])

sol = solve(prob,DFBDF(linsolve=LUFactorization()),maxiters=10000000,save_everystep=true,dtmax=1.0,verbose=false)
