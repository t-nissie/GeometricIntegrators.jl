GeometricIntegrators.jl
==============================
Geometric integrators (Euler, velocity-Verlet, position-Verlet and leapfrog) in Julia

You can use this module of GeometricIntegrators with Julia-0.4.x and Julia-0.5.x.

Its homepage is https://github.com/t-nissie/GeometricIntegrators.jl .

## Setup
There are no dependencies, just install with

    Pkg.clone("git://github.com/t-nissie/GeometricIntegrators.jl.git")

## Example
### In-line example
An example is in the end of `src/GeometricIntegrators.jl`.
If you are using Julia-0.5 or higher and
Winston (https://github.com/nolta/Winston.jl),

    $ cd src/
    $ julia GeometricIntegrators.jl

gives you an example plot of `GeometricIntegrators.eps`.

### Harmonic oscillator
`src/harmonic_oscillator.jl` is an example of one-dimensional harmonic oscillator.
Execute `harmonic_oscillator.jl`, then you will get `harmonic_oscillator_??.eps`.

### Two body problem
`src/twobody.jl` is an example of two body problem.
Execute `twobody.jl`, then you will get `twobody_??.eps` as shown in Fig. 1.

![twobody](https://raw.githubusercontent.com/t-nissie/GeometricIntegrators.jl/master/docs/twobody.jpg "two body problem")

FIG. 1. Trajectories of two bodies with masses M=2 and m=1.

## Usage
Prepare arrays of functions qdot and pdot.
