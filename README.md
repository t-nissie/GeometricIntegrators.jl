GeometricIntegrators.jl
=======================
[Geometric integrators](https://en.wikipedia.org/wiki/Geometric_integrator)
(Euler, symplectic-Euler, velocity-Verlet, position-Verlet and leapfrog) in Julia

You can use this module of GeometricIntegrators with Julia-1.4.

Its homepage is https://github.com/t-nissie/GeometricIntegrators.jl .

GeometricIntegrators.jl is under continuous integration at Travis CI:
[![Build Status](https://travis-ci.org/t-nissie/GeometricIntegrators.jl.svg?branch=master](https://travis-ci.org/t-nissie/GeometricIntegrators.jl)

## Setup
There are no dependencies, just install with

    using Pkg; Pkg.add(PackageSpec(url="https://github.com/t-nissie/GeometricIntegrators.jl", rev="master"));

## Example
### In-line example
An example is in the end of `src/GeometricIntegrators.jl`.
If you are using Julia-1.4 and
Winston (https://github.com/nolta/Winston.jl),

    $ cd src/
    $ julia GeometricIntegrators.jl

gives you an example plot of `GeometricIntegrators.eps`.

### Harmonic oscillator
`test/harmonic_oscillator.jl` is an example of one-dimensional harmonic oscillator.
Execute `harmonic_oscillator.jl`, then you will get `harmonic_oscillator_??.eps`.

### Two body problem
`test/twobody.jl` is an example of two body problem.
Execute `twobody.jl`, then you will get `twobody_??.eps` as shown in Fig. 1.

![twobody](https://raw.githubusercontent.com/t-nissie/GeometricIntegrators.jl/master/docs/twobody.jpg "two body problem")

FIG. 1. Trajectories of two bodies with masses M=2 and m=1.

### Test
Clone, test, then remove this package:

    using Pkg; Pkg.add(PackageSpec(url="https://github.com/t-nissie/GeometricIntegrators.jl", rev="master")); Pkg.test("GeometricIntegrators"); Pkg.rm("GeometricIntegrators")

## Usage
Prepare arrays of functions qdot and pdot.
