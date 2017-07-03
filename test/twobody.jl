#!/usr/bin/env julia
##
using GeometricIntegrators
using Base.Test

M = 2.0
m = 1.0
G = 1.0
qdot = [p ->  p[1]/M,
        p ->  p[2]/M,
        p ->  p[3]/m,
        p ->  p[4]/m]
pdot = [q -> G*M*m*(q[3]-q[1]) / ((q[3]-q[1])^2+(q[4]-q[2])^2)^1.5,
        q -> G*M*m*(q[4]-q[2]) / ((q[3]-q[1])^2+(q[4]-q[2])^2)^1.5,
        q -> G*M*m*(q[1]-q[3]) / ((q[3]-q[1])^2+(q[4]-q[2])^2)^1.5,
        q -> G*M*m*(q[2]-q[4]) / ((q[3]-q[1])^2+(q[4]-q[2])^2)^1.5]
q0 = [-3.0, 0.0,  6.0,0.0]
p0 = [ 0.0,-0.2,  0.0,0.2]
t0 = 0.0
h  = 0.1
n  = 400
p_half = p0 + h/2 .* map(f -> f(q0), pdot)

t_eu,q_eu,p_eu = time_evolution(          euler, qdot, pdot, q0, p0, t0, h, n)
t_vv,q_vv,p_vv = time_evolution(velocity_verlet, qdot, pdot, q0, p0, t0, h, n)
t_pv,q_pv,p_pv = time_evolution(position_verlet, qdot, pdot, q0, p0, t0, h, n)
t_lf,q_lf,p_lf = time_evolution(       leapfrog, qdot, pdot, q0, p_half, t0, h, n)

@test q_eu[1,  2] ≈ q_vv[1,  2] atol=1.0e-3
@test q_eu[2,  2] ≈ q_vv[2,  2] atol=1.0e-3
@test q_eu[3,  2] ≈ q_vv[3,  2] atol=1.0e-3
@test q_eu[4,  2] ≈ q_vv[4,  2] atol=1.0e-3

@test q_pv[1,n+1] ≈ q_vv[1,n+1] atol=1.0e-4
@test q_pv[2,n+1] ≈ q_vv[2,n+1] atol=1.0e-4
@test q_pv[3,n+1] ≈ q_vv[3,n+1] atol=1.0e-4
@test q_pv[4,n+1] ≈ q_vv[4,n+1] atol=1.0e-4

@test q_lf[1,n+1] ≈ q_vv[1,n+1] atol=1.0e-13
@test q_lf[2,n+1] ≈ q_vv[2,n+1] atol=1.0e-13
@test q_lf[3,n+1] ≈ q_vv[3,n+1] atol=1.0e-13
@test q_lf[4,n+1] ≈ q_vv[4,n+1] atol=1.0e-13

using Winston
plot(q_pv[1,:],q_pv[2,:],"r-", q_pv[1,1:20:n],q_pv[2,1:20:n],"ro",
     q_pv[3,:],q_pv[4,:],"g--",q_pv[3,1:20:n],q_pv[4,1:20:n],"gs")
xlim(-5.0,7.0)
ylim(-6.0,6.0)
savefig("twobody_pv.eps")
plot(q_lf[1,:],q_lf[2,:],"r-", q_lf[1,1:20:n],q_lf[2,1:20:n],"ro",
     q_lf[3,:],q_lf[4,:],"g--",q_lf[3,1:20:n],q_lf[4,1:20:n],"gs")
xlim(-5.0,7.0)
ylim(-6.0,6.0)
savefig("twobody_lf.eps")
