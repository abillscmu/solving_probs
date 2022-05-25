using OrdinaryDiffEq
using LinearSolve
using IterativeSolvers
using ComponentArrays
using Parameters

p = ComponentArray(a=0.04,b=1e4,c=3e7,d=1.0)
u0 = ComponentArray(x1=1.0,x2=0.0,x3=0.0)

function f(resid,du,u,p,t)
    @unpack a,b,c,d=p
    resid[1] = - a*u[1] + b*u[2]*u[3] - du[1]
    resid[2] = + a*u[1] - c*u[2]^2 - b*u[2]*u[3] - du[2]
    resid[3] = u[1] + u[2] + u[3] - d
end

du0 = similar(u0)
fill!(du0,0.0)

tspan = (0.0,1e5)


prob = DAEProblem(f,du0,u0,tspan,p,differential_vars=[true,true,false])

sol = solve(prob,DFBDF(linsolve=LUFactorization()))


#component arrays -- DFBDF solver
