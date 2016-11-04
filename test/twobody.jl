#!/usr/bin/env julia
##
using MultipleTimeStepIntegrators

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

t,q,p = time_evolution(position_verlet, qdot, pdot, q0, p0, t0, h, n)

using Winston
plot(q[1,:],q[2,:],"r-", q[1,1:20:n],q[2,1:20:n],"ro",
     q[3,:],q[4,:],"g--",q[3,1:20:n],q[4,1:20:n],"gs")
xlim(-5.0,7.0)
ylim(-6.0,6.0)
savefig("mts3.eps")
