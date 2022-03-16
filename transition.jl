# Motoaki Takahashi 

# Julia 1.7.2

# I compute a transition path from a given initial condition to a steady state 

using LinearAlgebra, Statistics, Compat, Plots
using CSV
using RData
using DataFrames
using SpecialFunctions
using Plots: length
using CSVFiles
using Query
using QuerySQLite
using RData
using Missings

# First I need to compute a steady steat toward which a dynamic path converges

# suppose the economy converges to a steady state in 100 periods from an initial period

T = 100

# the number of places 
N = 2

# the number of cohorts in a given period
C = 8

# the number of races 
R = 2

# for parameters, columns generically represent periods,
# rows are races/places 

# migration elesticity
ν = fill(2.0, 1, T) # this follows Caliendo, Opromolla, Parro, and Sforza

# elasticity of substitution between ages 
σ_0 = fill(5.0, 1, T) # this follows Ottaviano and Peri, Manacorda et al, and Yuta's JMP

# elasticity of substitution between races (within ages)
σ_1 = fill(11.0, 1, T) # this follows Boustan

# Cobb-Douglas share of housing
γ = fill(0.25, 1, T)

# amenities 
B_period = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
# the amenity for the youngest race 1 in place 1, ..., the amenity for the oldest race 1 in place 1, ...
# the amenity for the youngest race 1 in place 2, ..., the amenity for the oldest race 1 in place 2, ...
# the amenity for the youngest race 2 in place 1, ..., the amenity for the oldest race 2 in place 1, ...
# the amenity for the youngest race 2 in place 2, ..., the amenity for the oldest race 2 in place 2, ...

B = repeat(B_period, 1, T)

# survival probabilities
s_period = [       1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                   1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
# the survival probability for the youngest race 1, ..., the survival probability for the second oldest race 1, ...
# the survival probability for the youngest race 2, ..., the survival probability for the second oldest race 2, ...
s = repeat(s, 1, T)

# migration costs 
# note that the oldest can no longer migrate

τ_period = [     0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
                 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0,
                 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
                 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0]

# the 1st row:
# the youngest race 1's migration cost from place 1 to place 1, the youngest race 1's migration cost from place 1 to place 2, ...,
# the second oldest race 1's migration cost from place 1 to place 1, the second oldest race 1's migration cost from place 1 to place 2, ...,
# ...,
# the last row:
# the youngest race 2's migration cost from place 2 to place 1, the youngest race 2's migration cost from place 2 to place 2, ...,
# the second oldest race 2's migration cost from place 2 to place 1, the second oldest race 2's migration cost from place 2 to place 2

τ = repeat(τ_period, 1, T)


# place-specific shifter of rent 
r_bar_period = [0.1, 0.1]

r_bar = repeat(r_bar_period, 1, T)

# place-specific elasticity of rent 
η_period = [0.5, 0.5]

η = repeat(η, 1, T)

# productivity
A_period = [1.0, 1.0]
# the productivity in place 1, the productivity in place 2

A = repeat(A, 1, T)

# cohort-specific productivity 
κ_0_period = [1/7, 1/7, 1/7, 1/7, 1/7, 1/7, 1/7,
              1/7, 1/7, 1/7, 1/7, 1/7, 1/7, 1/7]
# the relative productivity of the second youngest in place 1, ..., the relative productivity of the oldest in place 1,
# the relative productivity of the second youngest in place 2, ..., the relative productivity of the oldest in place 2

κ_0 = repeat(κ_0_period, 1, T)

# the following is copied from steady_state.jl 

####################################################################

