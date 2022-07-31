GeometricIntegrators.jl
=======================
[Geometric integrators](https://en.wikipedia.org/wiki/Geometric_integrator)
(Euler, symplectic-Euler, velocity-Verlet, position-Verlet and leapfrog) in Julia

You can use this module of GeometricIntegrators with Julia-1.6 and Julia-1.7.

Its homepage is https://github.com/t-nissie/GeometricIntegrators.jl .

GeometricIntegrators.jl is under continuous integration at Travis CI:
[![Build Status](https://travis-ci.com/t-nissie/GeometricIntegrators.jl.svg?branch=master)](https://travis-ci.com/github/t-nissie/GeometricIntegrators.jl)

Coverrage report:
[![Coverage Status](https://coveralls.io/repos/github/t-nissie/GeometricIntegrators.jl/badge.svg?branch=master)](https://coveralls.io/github/t-nissie/GeometricIntegrators.jl?branch=master)

## Setup
There are no dependencies, just install with

    julia> using Pkg; Pkg.add(PackageSpec(url="https://github.com/t-nissie/GeometricIntegrators.jl", rev="master"));

You may also need to install Winston package.

    julia> using Pkg; Pkg.add("Winston")

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

    $ cd test/
    $ julia --project harmonic_oscillator.jl

### Two body problem
`test/twobody.jl` is an example of two body problem.
Execute `twobody.jl`, then you will get `twobody_??.eps` as shown in Fig. 1.

     $ cd test/
     $ julia --project twobody.jl

![twobody](https://raw.githubusercontent.com/t-nissie/GeometricIntegrators.jl/master/docs/twobody.jpg "two body problem")

FIG. 1. Trajectories of two bodies with masses M=2 and m=1.

## Test
### Test as a user
Clone, test, then remove this package:

    julia> using Pkg; Pkg.add(PackageSpec(url="https://github.com/t-nissie/GeometricIntegrators.jl", rev="master")); Pkg.test("GeometricIntegrators"); Pkg.rm("GeometricIntegrators")


### Test as a developer
From Pkg REPL-mode:

    $ julia --project
    julia> ]
    (GeometricIntegrators) pkg> test

From command line:

    $ julia --project -e 'using Pkg; Pkg.test()'


## Usage
Prepare arrays of functions qdot and pdot.
