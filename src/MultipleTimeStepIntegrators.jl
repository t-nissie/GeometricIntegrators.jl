#
# MultipleTimeStepIntegrators.jl: multiple-time-step (MTS) integrators in Julia
# Copyright (C) 2016 Takeshi Nishimatsu
#
# License: GPLv3
# Setup: Pkg.clone("git://github.com/t-nissie/MultipleTimeStepIntegrators.jl.git")
# Example: An example is in the end of this file.
#          If you are using Julia-0.5 or higher and
#          Winston (https://github.com/nolta/Winston.jl),
#          "julia MultipleTimeStepIntegrators.jl" gives you
#          a plot of MultipleTimeStepIntegrators.eps .
# References:
# Ref1. 奥村久士:『分子動力学シミュレーションにおける温度・圧力制御
#       第2回:シンプレクティック分子動力学法と能勢・ポアンカレ熱浴』，
#       分子シミュレーション研究会会誌《アンサンブル》
#       Vol.11, No.1, January 2009（通巻45号）(in Japanese).
# Ref2. M. Tuckerman, B. J. Berne, and G. J. Martyna:
#       "Reversible multiple time scale molecular dynamics",
#       J. Chem. Phys. 97, 1990 (1992).
# Ref3. Paul F. Batcho and Tamar Schlick:
#       "Special stability advantages of position-Verlet over
#       velocity-Verlet in multiple-time step integration",
#       J. Chem. Phys. 115, 4019 (2001).
# Ref4. https://en.wikipedia.org/wiki/Leapfrog_integration
# Ref5. https://github.com/timothyrenner/RungeKutta.jl
##
module MultipleTimeStepIntegrators

function           euler{T<:AbstractFloat}(qdot,   # It was qdot::Array{Function,1}.
                                           pdot,
                                              q::Array{T,1},
                                              p::Array{T,1},
                                              t::AbstractFloat,
                                              h::AbstractFloat)
    q_next = q      + h   .* map(f -> f(p), qdot)
    p_next = p      + h   .* map(f -> f(q), pdot)
    return t + h, q_next, p_next
end

function velocity_verlet{T<:AbstractFloat}(qdot,
                                           pdot,
                                              q::Array{T,1},
                                              p::Array{T,1},
                                              t::AbstractFloat,
                                              h::AbstractFloat)
    p_half = p      + h/2 .* map(f -> f(q),      pdot)
    q_next = q      + h   .* map(f -> f(p_half), qdot)
    p_next = p_half + h/2 .* map(f -> f(q_next), pdot)
    return t + h, q_next, p_next
end

function position_verlet{T<:AbstractFloat}(qdot,
                                           pdot,
                                              q::Array{T,1},
                                              p::Array{T,1},
                                              t::AbstractFloat,
                                              h::AbstractFloat)
    q_half = q      + h/2 .* map(f -> f(p),      qdot)
    p_next = p      + h   .* map(f -> f(q_half), pdot)
    q_next = q_half + h/2 .* map(f -> f(p_next), qdot)
    return t + h, q_next, p_next
end

function        leapfrog{T<:AbstractFloat}(qdot,
                                           pdot,
                                              q::Array{T,1},
                                              p::Array{T,1},
                                              t::AbstractFloat,
                                              h::AbstractFloat)
    q_next = q + h .* map(f -> f(p), qdot) + h^2 .* map(f -> f(q), pdot) / 2
    p_next = p + h .* (map(f -> f(q), pdot) + map(f -> f(q_next), pdot)) / 2
    return t + h, q_next, p_next
end

function time_evolution{T<:AbstractFloat}(method,
                                          qdot,
                                          pdot,
                                            q0::Array{T,1},
                                            p0::Array{T,1},
                                            t0::AbstractFloat,
                                             h::AbstractFloat,
                                             n::Integer)
    q = zeros(length(q0), n+1)
    p = zeros(length(p0), n+1)
    t = zeros(1, n+1)
    q[:,1] = q0
    p[:,1] = p0
    t[1] = t0
    for ii in 2:(n+1)
        t[ii], q[:,ii], p[:,ii] = method(qdot, pdot, q[:,ii-1], p[:,ii-1], t[ii-1], h)
    end
    return t,q,p
end

export euler, velocity_verlet, position_verlet, leapfrog, time_evolution
end

if isdefined(:PROGRAM_FILE) && PROGRAM_FILE == basename(@__FILE__)
   using MultipleTimeStepIntegrators

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
    qdot = [p ->  p[1]/m,
            p ->  p[2]/m]
    pdot = [q -> -q[1]*k,
            q -> -q[2]*k]
    q0 = [0.0,2.0]
    p0 = [1.0,1.0]
    t0 = 0.0
    h  = 0.01
    n  = 1000

    t,q,p = time_evolution(position_verlet, qdot, pdot, q0, p0, t0, h, n)
    e = energy(q,p,m,k,n)

    using Winston   # It takes some time.
    plot(t, q[1,:], t, p[1,:], t, q[2,:], t, p[2,:], t, e[:])
    #ylim(-2.5,2.5)
    savefig("MultipleTimeStepIntegrators.eps")
end
