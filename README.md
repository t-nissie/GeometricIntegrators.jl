MultipleTimeStepIntegrators.jl
==============================
multiple time step (MTS) integrators (Euler, velocity-Verlet, position-Verlet and leapfrog) in Julia

You can use this module of MultipleTimeStepIntegrators with Julia-0.4.x and Julia-0.5.x.

Its homepage is https://github.com/t-nissie/MultipleTimeStepIntegrators.jl .

## Setup
There are no dependencies, just install with

    Pkg.clone("git://github.com/t-nissie/MultipleTimeStepIntegrators.jl.git")

## Example
### In-line example
An example is in the end of `src/MultipleTimeStepIntegrators.jl`.
If you are using Julia-0.5 or higher and
Winston (https://github.com/nolta/Winston.jl),

    $ cd src/
    $ julia MultipleTimeStepIntegrators.jl

gives you an example plot of `MultipleTimeStepIntegrators.eps`.

### Two body problem
`src/twobody.jl` is an example of two body problem.
Execute `twobody.jl`, then you will get `twobody.eps` as shown in Fig. 1.

![twobody](https://raw.githubusercontent.com/t-nissie/MultipleTimeStepIntegrators.jl/master/docs/twobody.jpg "two body problem")

FIG. 1. Trajectories of two bodies with masses M=2 and m=1.

## Usage
Prepare arrays of functions qdot and pdot.
