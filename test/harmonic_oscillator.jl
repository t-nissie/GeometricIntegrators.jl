#!/usr/bin/env julia
##
using GeometricIntegrators
using Base.Test

function energy{T<:AbstractFloat}(q::Array{T,2},
                                  p::Array{T,2}, m, k, n::Integer)
    e = zeros(Float64, n+1)
    for i in 1:n+1
        e[i] = dot(p[:,i],p[:,i])/2/m +
               dot(q[:,i],q[:,i])*k/2
    end
    return e
end

# H = p^2/(2m) + k*q^2/2
m = 1.0
k = 1.0
qdot = [p ->  p[1]/m]
pdot = [q -> -q[1]*k]
q0 = [0.0]
p0 = [2.0]
t0 = 0.0
h  = 0.01
n  = 1000
p_half = p0 + h/2 .* map(f -> f(q0), pdot)

t_eu,q_eu,p_eu = time_evolution(          euler, qdot, pdot, q0, p0,     t0, h, n)
energy_eu = energy(q_eu,p_eu,m,k,n)
t_vv,q_vv,p_vv = time_evolution(velocity_verlet, qdot, pdot, q0, p0,     t0, h, n)
energy_vv = energy(q_vv,p_vv,m,k,n)
t_pv,q_pv,p_pv = time_evolution(position_verlet, qdot, pdot, q0, p0,     t0, h, n)
energy_pv = energy(q_pv,p_pv,m,k,n)
t_lf,q_lf,p_lf = time_evolution(       leapfrog, qdot, pdot, q0, p_half, t0, h, n)

@test q_eu[1,n+1] ≈ q_vv[1,n+1] atol=1.0e-1
@test q_pv[1,n+1] ≈ q_vv[1,n+1] atol=1.0e-4
@test q_lf[1,n+1] ≈ q_vv[1,n+1] atol=1.0e-14

using Winston
plot(t_eu, q_eu[1,:], t_eu, p_eu[1,:], t_eu, energy_eu[:])
ylim(-2.5,2.5)
savefig("harmonic_oscillator_eu.eps")
plot(t_vv, q_vv[1,:], t_vv, p_vv[1,:], t_vv, energy_vv[:])
ylim(-2.5,2.5)
savefig("harmonic_oscillator_vv.eps")
plot(t_pv, q_pv[1,:], t_pv, p_pv[1,:], t_pv, energy_pv[:])
ylim(-2.5,2.5)
savefig("harmonic_oscillator_pv.eps")
plot(t_lf, q_lf[1,:], t_lf, p_lf[1,:])
ylim(-2.5,2.5)
savefig("harmonic_oscillator_lf.eps")
