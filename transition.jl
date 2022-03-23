# Motoaki Takahashi 

# Julia 1.7.2

# March 2022

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

# suppose the economy converges to a steady state in 100 periods from an initial period
# In period 0, the economy is at the steady state.
# In period 1, a shock hits the economy.
# From period 2 to period 99, the economy is on a transition path to the steady state.
# In period 100, the economy reaches the steady state.

T = 101

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
s = repeat(s_period, 1, T)

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

η = repeat(η_period, 1, T)

# productivity
A_period = [4.5, 5.5]
# the productivity in place 1, the productivity in place 2

A = repeat(A_period, 1, T)

# cohort-specific productivity 
κ_0_period = [1/7, 1/7, 1/7, 1/7, 1/7, 1/7, 1/7,
              1/7, 1/7, 1/7, 1/7, 1/7, 1/7, 1/7]
# the relative productivity of the second youngest in place 1, ..., the relative productivity of the oldest in place 1,
# the relative productivity of the second youngest in place 2, ..., the relative productivity of the oldest in place 2

κ_0 = repeat(κ_0_period, 1, T)

# race-specific productivity (within cohorts)
κ_1_period = [ 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2,
               1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2,
               1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2,
               1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2]
# the relative productivity of race1 within the oldest in place 1, ..., the relative productivity of race1 within the second youngest in place 1,
# the relative productivity of race1 within the oldest in place 2, ..., the relative productivity of race1 within the second youngest in place 2,
# the relative productivity of race2 within the oldest in place 1, ..., the relative productivity of race2 within the second youngest in place 1,
# the relative productivity of race2 within the oldest in place 2, ..., the relative productivity of race2 within the second youngest in place 2

κ_1 = repeat(κ_1_period, 1, T)


# fertility per cohort-race-place 
# I guess I need to restrict the value of α to get a steady state, 
α_period = [    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,]
# the fertility of the second youngest race 1 in place 1, ..., the fertility of the oldest race 1 in place 1,
# the fertility of the second youngest race 2 in place 1, ..., the fertility of the oldest race 2 in place 1,
# the fertility of the second youngest race 2 in place 1, ..., the fertility of the oldest race 2 in place 1,
# the fertility of the second youngest race 2 in place 2, ..., the fertility of the oldest race 2 in place 2

α = repeat(α_period, 1, T)

# immgrants from abroad 
M_period = [    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
# the immigrants of the youngest race 1 in place 1, ..., the immigrants of the oldest race 1 in place 1,
# the immigrants of the youngest race 1 in place 2, ..., the immigrants of the oldest race 1 in place 2,
# ...,
# the immigrants of the youngest race 2 in place 2, ..., the immigrants of the oldest race 2 in place 2.

M = repeat(M_period, 1, T)


# the following functions are the same as in steady_state_from_young_to_old.jl 

function utility(w, r, B, R, N, C, γ)
        temp0 = kron(r .^ γ, fill(1.0, (C-1)))
        temp1 = repeat(temp0, R)
        
        
        temp2 = log.(w./ temp1)
        u = []
        
        for i in 1:(R*N)
                u = [u; 0.0; temp2[((i-1)*(C-1)+1):(i*(C-1))] ]
        end
        
        u = u + log.(B)
        
        return u
end

function expected_value(u, s, ν, τ, R, N, C)
    # transform vectors to a matrices
    τ_mat = reshape(τ, (N*(C-1)), (R*N))
    τ_mat = τ_mat'

    u_mat = reshape(u, C, (R*N))
    u_mat = u_mat'

    s_mat = reshape(s, (C-1), R)
    s_mat = s_mat'

    V_mat = zeros((R*N), C)

    V_mat[:, C] =  u_mat[:, C]
    
    for i in 1:R
            # we fill in expected values from the second oldest cohort to the youngest cohort
            for k in 1:(C-1)
                    for j in 1:N                  
                            V_mat[((i-1)*N+j), (C-k)] = (u_mat[((i-1)*N+j), (C-k)] 
                            + ν * log(sum(exp.(s_mat[i, (C-k)] * V_mat[((i-1)*N+1):(i*N),(C-k+1)] - τ_mat[((i-1)*N+j), ((C-k-1)*N+1):((C-k)*N)]) .^ (1/ν) ) ) )
                    end
            end
    end

    # return a vector, not a matrix
    output = reshape(V_mat', (R*N*C) )
    return output

end

function migration_rate(s, V, ν, τ, R, N, C)
    V_mat = reshape(V, C, (R*N))'
    s_mat = reshape(s, (C-1), R)'
    τ_mat = reshape(τ, (N*(C-1)), (R*N))'
    
    # first I make a (R*N*N) × (C-1) migration share matrix
    μ_mat = zeros((R*N*N), (C-1))
    
    for i in 1:R
            for k in 1:(C-1)
                    for j in 1:N
                    μ_mat[((i-1)*N*N + (j-1)*N + 1):((i-1)*N*N + j*N), k] = (exp.(s_mat[i, k] * V_mat[((i-1)*N+1):(i*N), (k+1)] - τ_mat[((i-1)*N + j), ((k-1)*N+1):(k*N)]) .^ (1/ν) 
                    ./ sum(exp.(s_mat[i, k] * V_mat[((i-1)*N+1):(i*N), (k+1)] - τ_mat[((i-1)*N + j), ((k-1)*N+1):(k*N)]) .^ (1/ν) ) )
                    end
            end
    
    end
    # μ_mat is:
    # the migration rate of the youngest race 1 from place 1 to place 1, ..., the migration rate of the second oldest race 1 from place 1 to place 1,
    # the migration rate of the youngest race 1 from place 1 to place 2, ..., the migration rate of the second oldest race 1 from place 1 to place 2,
    # ...,
    # the migration rate of the youngest race 2 from place 2 to place 1, ..., the migration rate of the second oldest race 2 from place 2 to place 1,
    # the migration rate of the youngest race 2 from place 2 to place 2, ..., the migration rate of the second oldest race 2 from place 2 to place 2.


    output = reshape(μ_mat', (R*N*N*(C-1)) )
    return output
end

function population(μ, s, L, α, C, R, N)
    μ_mat = reshape(μ, (C-1), (R*N*N))'
    s_mat = reshape(s, (C-1), R)'
    α_mat = reshape(α, C, (R*N))'
    L_mat_input = reshape(L, C, (R*N))'
    
    L_mat_output = zeros((R*N), C)
    
    # the number of the new-born
    L_mat_output[:, 1] = sum(α_mat .* L_mat_input, dims = 2)
    
    temp = kron(L_mat_input[:,1:(C-1)], ones(N)) .* kron(s_mat, ones(N*N)) .* μ_mat 
    
    for i in 1:R
            for k in 1:(C-1)
                    for j in 1:N
                            for l in 1:N
                                    L_mat_output[((i-1)*N+j), (k+1)] = L_mat_output[((i-1)*N+j), (k+1)] + temp[((i-1)*N*N + j + (l-1)*N) , k]
                            end
                    end
            end
    end
    
    output = reshape(L_mat_output', (R*N*C))

end

function agg_labor_cohort(L, σ_1, κ_1, N, R, C)
    κ_1_mat = reshape(κ_1, (C-1), (R*N))'
    L_mat = reshape(L, C, (R*N))'
    
    
    temp = (κ_1_mat .^ (1/σ_1)) .* (L_mat[:, 2:C] .^ ((σ_1 - 1)/σ_1))
    
    # note that the youngest cannot work
    L_cohort_mat = zeros(N, (C-1))
    
    for i in 1:N
            for j in 1:R
            L_cohort_mat[i,:] = L_cohort_mat[i,:] + temp[(i+(j-1)*N), :]
            end
    end
    
    L_cohort_mat = L_cohort_mat .^ (σ_1/(σ_1-1))
    
    output = reshape(L_cohort_mat', (C-1)*N)

    return output

end

function agg_labor_place(L_cohort, σ_0, κ_0, N, R, C)
    L_cohort_mat = reshape(L_cohort, (C-1), N)'
    κ_0_mat = reshape(κ_0, (C-1), N)'
    
    temp = κ_0_mat .^ (1/σ_0) .* L_cohort_mat .^ ((σ_0-1)/σ_0)
    L_place_mat = sum(temp, dims = 2) .^ (σ_0/(σ_0-1))

    return L_place_mat

end

function wage(σ_0, σ_1, A, κ_0, κ_1, N, R, C, L_place, L_cohort, L)
    L_cohort_mat = reshape(L_cohort, (C-1), N)'

    κ_0_mat = reshape(κ_0, (C-1), N)'
    κ_1_mat = reshape(κ_1, (C-1), (R*N))'
    
    L_mat = reshape(L, C, (R*N))'
    
    # wage_mat is of (N*R) × (C-1)
    wage_mat = (repeat(A, R, (C-1)) .* repeat(L_place, R, (C-1)) .^ (1/(σ_0-1)) .* repeat(κ_0_mat, R, 1) .^ (1/σ_0)
            .* repeat(L_cohort_mat, R, 1) .^ (-1/σ_0 + 1/(σ_1-1)) .* κ_1_mat .^ (1/σ_1) .* L_mat[:, 2:C] .^ (-1/σ_1))


    output = reshape(wage_mat', (R*N*(C-1)))

end

function rent(w, L, r_bar, η, γ, N, R, C)
    w_mat = reshape(w, (C-1), (N*R))'
    L_mat = reshape(L, C, (N*R))'
    temp = sum(w_mat .* L_mat[:, 1:(C-1)], dims = 2)
    output = zeros(N)
    
    for i in 1:N
            for j in 1:R 
                    output[i] = output[i] + temp[i+(j-1)*N]
            end
    end

    output = r_bar .* (γ * output) .^ η

    return output
end


# First I need to compute a steady steat toward which a dynamic path converges

# an initial guess for population over time 
L_in = fill(10, (R*N*C, T))
# rows are race-place-age-in-the-current-periods, columns are periods
# in each column (period), 
# the 1st row is the population of the youngest race 1 in place 1,
# the 2nd row is the population of the second youngest race 1 in place 1,
# ...,
# the 32nd row is the population of the oldest race 2 in place 2

tol = 10 ^ (-7)
maxit = 1000

# the dumpening parameter for iteration
λ = 0.5

function steady_state(tol, maxit, N, C, R, ν, σ_0, σ_1, γ, B, s, τ, r_bar, η, A, κ_0, κ_1, α, L_in, λ)
        L = L_in 
        V = []

        dif = 1.0
        count = 0

        while dif > tol && count < maxit
                L_cohort = agg_labor_cohort(L, σ_1, κ_1, N, R, C)
                L_place = agg_labor_place(L_cohort, σ_0, κ_0, N, R, C)

                w = wage(σ_0, σ_1, A, κ_0, κ_1, N, R, C, L_place, L_cohort, L)
                r = rent(w, L, r_bar, η, γ, N, R, C)
                u = utility(w, r, B, R, N, C, γ)
                V = expected_value(u, s, ν, τ, R, N, C)
                μ = migration_rate(s, V, ν, τ, R, N, C)
                L_new = population(μ, s, L, α, C, R, N)

                dif = maximum(abs.((L - L_new) ./ L))

                count = count + 1

                L = (1-λ) * L + λ * L_new
        
        end

        if dif < tol
                output = [L V]
        else 
                output = fill("not converged", length([L V]))
        end

        return output
end

LV_100 = steady_state(tol, maxit, N, C, R, ν[1, T], σ_0[1, T], σ_1[1, T], γ[1, T], B[:, T], s[:, T], τ[:, T], r_bar[:, T], η[:, T], A[:, T], κ_0[:, T], κ_1[:, T], α[:, T], L_in[:, T], λ)

L_mat_100 = reshape(LV_100[:, 1], C, R*N)'
V_mat_100 = reshape(LV_100[:, 2], C, R*N)'


maximum(abs.(L_mat_100[1, :] - L_mat_100[3, :]))


# Next I compute population backward from the last period T

L = fill(10.0, (R*N*C), T)

V = fill(0.0, (R*N*C), T)

# the last period is the steady state
L[:, T] = LV_100[:, 1]
V[:, T] = LV_100[:, 2]

# the initial period (0) is the steady state
L[:, 1] = LV_100[:, 1]
V[:, 1] = LV_100[:, 2]
# the "perceived" expected value before an unanticipated (MIT) shock happens in period 1

# let's compute the populations in period T-1
μ_T_minus_1 = migration_rate(s[:, T-1], V[:, T], ν[:, T-1], τ[:, T-1], R, N, C)

μ_mat_T_minus_1 = reshape(μ_T_minus_1, (C-1), (R*N*N))'

μ_mat_T_minus_1_wide = zeros(N*R, N*(C-1))


for i in 1:R
        for j in 1:(C-1)
                for k in 1:N
                        μ_mat_T_minus_1_wide[((i-1)*N+1):(i*N) , (j-1)*N+k] = μ_mat_T_minus_1[((i-1)*N*N+(k-1)*N+1):((i-1)*N*N+k*N), j]
                end

        end
end

s_T_minus_1 = s[:, T-1]
s_mat_T_minus_1 = reshape(s_T_minus_1, (C-1), R)'
s_mat_T_minus_1_wide = kron(s_mat_T_minus_1, ones(N, N))

s_μ_T_minus_1_wide = s_mat_T_minus_1_wide .* μ_mat_T_minus_1_wide

L_T_minus_1_mat = reshape(L[:, T-1], C, R*N)'
# zeros(R*N, C)

M_100_mat = reshape(M[:, T], C, R*N)'

# first I compute the populations of the youngest to the second oldest

for i in 1:R
        for j in 1:(C-1)
                # cohort c in period t-1 is t-1-c years old, cohort c in period 
                L_T_minus_1_mat[((i-1)*N+1):(i*N), j] = ( inv(s_μ_T_minus_1_wide[((i-1)*N+1):(i*N), ((j-1)*N+1):(j*N)])
                 * (L_mat_100[((i-1)*N+1):(i*N), (j+1)] - M_100_mat[((i-1)*N+1):(i*N), (j+1)]) )
                

        end

end

