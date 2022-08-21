# Motoaki Takahashi 

# Julia 1.7.2

# March, April, August 2022

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

# for the main desktop
cd("D:/onedrive/OneDrive - The Pennsylvania State University/dynamic/code")

# for the laptop
# cd("C:/Users/takah/OneDrive - The Pennsylvania State University/dynamic/code")

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
A_period = [6.0, 4.0]

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
# the fertility of the second youngest race 1 in place 2, ..., the fertility of the oldest race 1 in place 2, ???
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

    # the number of the new-borns
    L_mat_output[:, 1] = sum(α_mat .* L_mat_output, dims = 2)
    
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
    wage_mat = (repeat(A, R, (C-1)) .* repeat(L_place, R, (C-1)) .^ (1/(σ_0)) .* repeat(κ_0_mat, R, 1) .^ (1/σ_0)
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

function expected_value_transition(u, s, ν, τ, R, N, C, V_t_plus_1)
        # transform vectors to a matrices
        τ_mat = reshape(τ, (N*(C-1)), (R*N))
        τ_mat = τ_mat'
    
        u_mat = reshape(u, C, (R*N))
        u_mat = u_mat'
    
        s_mat = reshape(s, (C-1), R)
        s_mat = s_mat'
    
        V_mat_t_plus_1 = reshape(V_t_plus_1, C, (R*N))'
        V_mat = zeros(R*N, C)
    
        V_mat[:, C] =  u_mat[:, C]
        
        for i in 1:R
                # we fill in expected values from the second oldest cohort to the youngest cohort
                for k in 1:(C-1)
                        for j in 1:N                  
                                V_mat[((i-1)*N+j), (C-k)] = (u_mat[((i-1)*N+j), (C-k)] 
                                + ν * log(sum(exp.(s_mat[i, (C-k)] * V_mat_t_plus_1[((i-1)*N+1):(i*N),(C-k+1)] - τ_mat[((i-1)*N+j), ((C-k-1)*N+1):((C-k)*N)]) .^ (1/ν) ) ) )
                        end
                end
        end
    
        # return a vector, not a matrix
        output = reshape(V_mat', (R*N*C) )
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

tol = 10 ^ (-8)
maxit = 1000

# the dumpening parameter for iteration
λ = 0.5

function steady_state(tol, maxit, N, C, R, ν, σ_0, σ_1, γ, B, s, τ, r_bar, η, A, κ_0, κ_1, α, L_in, λ)
        L = zeros(R*N*C)
        L[:, :] = L_in[:, :] 
        V = []
        w = []
        r = []
        real_wage = []

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

        real_wage = w ./ (repeat(kron(r, ones((C-1))), R) .^ γ)
        real_wage = [real_wage; zeros(R*N)]

        if dif < tol
                output = [L V real_wage]
        else 
                output = fill("not converged", length([L V real_wage]))
        end

        return output
end

# compare steady states for various survival probabilities 

s_per_range = 0.5:0.01:1.0

K = length(s_per_range)

V_cs_s = zeros(N*R*C, K)
L_cs_s = zeros(N*R*C, K)
real_wage_cs_s = zeros(N*R*(C-1), K)
μ_cs_s = zeros(R*N*N*(C-1), K)

for i in 1:K
        s_per = [1.0, 1.0, 1.0, s_per_range[i], s_per_range[i], s_per_range[i], s_per_range[i],
                 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

        temp = steady_state(tol, maxit, N, C, R, ν[1, T], σ_0[1, T], σ_1[1, T], γ[1, T], B[:, T], s_per, τ[:, T], r_bar[:, T], η[:, T], A[:, T], κ_0[:, T], κ_1[:, T], α[:, T], L_in[:, T], λ)

        L_cs_s[:, i] = temp[:, 1]
        V_cs_s[:, i] = temp[:, 2]
        μ_cs_s[:, i] = migration_rate(s_per, temp[:, 2], ν[1, T], τ[:, T], R, N, C)
        
        real_wage_cs_s[:, i] = temp[1:(N*R*(C-1)), 3]
end

# expected value
# plot(s_per_range, V_cs_s[1, :], label = "age 0 of race 1 in location 1", legend = :bottomright)
# plot!(s_per_range, V_cs_s[C+1, :], label = "age 0 of race 1 in location 2")
# plot!(s_per_range, V_cs_s[2*C+1, :], label = "age 0 of race 2 in location 1")
# plot!(s_per_range, V_cs_s[3*C+1, :], label = "age 0 of race 2 in location 2")
# xlabel!("survival probability of race 1")
# ylabel!("expected value")
# savefig("../output/steady_state/suv_prob_exp_value.pdf")

# population 
# plot(s_per_range, L_cs_s[1, :], label = "age 0 of race 1 in location 1", legend = :bottomright)
# plot!(s_per_range, L_cs_s[2*C+1, :], label = "age 0 of race 2 in location 1")
# xlabel!("survival probability of race 1")
# ylabel!("population")
# savefig("../output/steady_state/suv_prob_pop_plc1.pdf")

# migration rate 
# plot(s_per_range, μ_cs_s[3, :], label = "age 2 of race 1 from 1 to 1", legend = :bottomright)
# plot!(s_per_range, μ_cs_s[4*(C-1)+3, :], label = "age 2 of race 2 from 1 to 1", legend = :bottomright)
# xlabel!("survival probability of race 1")
# ylabel!("migration rate")
# savefig("../output/steady_state/suv_prob_mig_1to1.pdf")


# plot(s_per_range, μ_cs_s[2*(C-1)+3, :], label = "age 3 of race 1 from 2 to 1", legend = :bottomright)
# plot!(s_per_range, μ_cs_s[4*(C-1)+2*(C-1)+3, :], label = "age 3 of race 2 from 2 to 1")
# xlabel!("survival probability of race 1")
# ylabel!("migration rate")
# savefig("../output/steady_state/suv_prob_mig_2to1.pdf")

# prepare for comparative transition paths

LV_100 = steady_state(tol, maxit, N, C, R, ν[1, T], σ_0[1, T], σ_1[1, T], γ[1, T], B[:, T], s[:, T], τ[:, T], r_bar[:, T], η[:, T], A[:, T], κ_0[:, T], κ_1[:, T], α[:, T], L_in[:, T], λ)

V_mat_100_2 = reshape(LV_100[:, 2], C, R*N)'
L_mat_100_2 = reshape(LV_100[:, 1], C, R*N)'

L_mat_100 = reshape(LV_100[:, 1], C, R*N)'
V_mat_100 = reshape(LV_100[:, 2], C, R*N)'
real_wage_mat_100 = reshape(LV_100[1:N*R*(C-1), 3], C-1, R*N)'

maximum(abs.(L_mat_100[1, :] - L_mat_100[3, :]))

# I compute populations forward, given expected values.
# Then given populations, I update expected values.

λ_2 = 1.0

L_in = fill(10.0, (R*N*C), T)

V_in = fill(1.0, (R*N*C), T)
V_in_2 = fill(2.0, (R*N*C), T)

# the last period is the steady state
L_in[:, T] = LV_100[:, 1]
V_in[:, T] = LV_100[:, 2]
V_in_2[:, T] = LV_100[:, 2]

# the initial period (0) is the steady state
L_in[:, 1] = LV_100[:, 1]
V_in[:, 1] = LV_100[:, 2]
V_in_2[:, 1] = LV_100[:, 2]
# the "perceived" expected value before an unanticipated (MIT) shock happens in period 1


function transition_path(R, C, N, T, λ_2, L_in, V_in, V_in_2, s, ν, τ, M, σ_0, σ_1, κ_0, κ_1, A, r_bar, η, γ, B, tol, maxit)
        L = L_in
        V = V_in 
        V_new = V_in_2
        dif = 1.0
        count = 0
        w = zeros(R*N*(C-1), T-1)
        r = zeros(N, T-1)


        while dif > tol && count < maxit

                for i in 1:T-2
                        if i == 1

                                                        # given a value matrix, compute populations forward
                        L_t = L[:, i]
                        s_t = s[:, i]
                        ν_t = ν[1, i]
                        τ_t = τ[:, i]
                        V_t_plus_1 = V[:, i+1]
                        α_t_plus_1 = α[:, i+1]
                        M_t_plus_1 = M[:, i+1]
                        μ_t = migration_rate(s_t, V_t_plus_1, ν_t, τ_t, R, N, C)
                        L_t_plus_1 = population(μ_t, s_t, L_t, α_t_plus_1, C, R, N) + M_t_plus_1
                
                        # given a population matrix, compute wages and rents. Then compute period utilities.
                        # Then update expected values.
                
                        σ_1_t_plus_1 = σ_1[1, i+1]
                        κ_1_t_plus_1 = κ_1[:, i+1]
                        σ_1_t = σ_1[1, i]
                        κ_1_t = κ_1[:, i]
                
                        L_cohort_t_plus_1 = agg_labor_cohort(L_t_plus_1, σ_1_t_plus_1, κ_1_t_plus_1, N, R, C)
                        L_cohort_t = agg_labor_cohort(L_t, σ_1_t, κ_1_t, N, R, C)
                
                        σ_0_t_plus_1 = σ_0[1, i+1]
                        κ_0_t_plus_1 = κ_0[:, i+1]
                        σ_0_t = σ_0[1, i]
                        κ_0_t = κ_0[:, i]
                
                
                        L_place_t_plus_1 = agg_labor_place(L_cohort_t_plus_1, σ_0_t_plus_1, κ_0_t_plus_1, N, R, C)
                        L_place_t = agg_labor_place(L_cohort_t, σ_0_t, κ_0_t, N, R, C)
                
                        A_t_plus_1 = A[:, i+1]
                        A_t = A[:, i]
                
                        w_t_plus_1 = wage(σ_0_t_plus_1, σ_1_t_plus_1, A_t_plus_1, κ_0_t_plus_1, κ_1_t_plus_1, N, R, C, L_place_t_plus_1, L_cohort_t_plus_1, L_t_plus_1)
                        w_t = wage(σ_0_t, σ_1_t, A_t, κ_0_t, κ_1_t, N, R, C, L_place_t, L_cohort_t, L_t)
                
                        r_bar_t_plus_1 = r_bar[:, i+1]
                        r_bar_t = r_bar[:, i]
                        γ_t_plus_1 = γ[1, i+1]
                        γ_t = γ[1, i]
                        η_t_plus_1 = η[1, i+1]
                        η_t = η[1, i]
                
                
                
                        r_t_plus_1 = rent(w_t_plus_1, L_t_plus_1, r_bar_t_plus_1, η_t_plus_1, γ_t_plus_1, N, R, C)
                        r_t = rent(w_t, L_t, r_bar_t, η_t, γ_t, N, R, C)
                
                        B_t_plus_1 = B[:, i+1]
                        B_t = B[:, i]
                
                        u_t_plus_1 = utility(w_t_plus_1, r_t_plus_1, B_t_plus_1, R, N, C, γ_t_plus_1)
                        u_t = utility(w_t, r_t, B_t, R, N, C, γ_t)
                
                        s_t_plus_1 = s[:, i+1]
                        τ_t_plus_1 = τ[:, i+1]
                        V_t_plus_2 = V[:, i+2]
                        ν_t_plus_1 = ν[1, i+1]

                        s_t = s[:, i]
                        τ_t = τ[:, i]
                        ν_t = ν[1, i]

                        V_t_plus_1 = expected_value_transition(u_t_plus_1, s_t_plus_1, ν_t_plus_1, τ_t_plus_1, R, N, C, V_t_plus_2)
                        V_t = expected_value_transition(u_t, s_t, ν_t, τ_t, R, N, C, V_t_plus_1)
                
                        V_new[:, i+1] = V_t_plus_1
                        V_new[:, i] = V_t
                        
                        L[:, i+1] = L_t_plus_1  
                        w[:, i] = w_t
                        w[:, i+1] = w_t_plus_1
                        r[:, i] = r_t 
                        r[:, i+1] = r_t_plus_1 

                        


                
    
                        
                        else

                        # given a value matrix, compute populations forward
                        L_t = L[:, i]
                        s_t = s[:, i]
                        ν_t = ν[1, i]
                        τ_t = τ[:, i]
                        V_t_plus_1 = V[:, i+1]
                        α_t_plus_1 = α[:, i+1]
                        M_t_plus_1 = M[:, i+1]
                        μ_t = migration_rate(s_t, V_t_plus_1, ν_t, τ_t, R, N, C)
                        L_t_plus_1 = population(μ_t, s_t, L_t, α_t_plus_1, C, R, N) + M_t_plus_1
                
                        # given a population matrix, compute wages and rents. Then compute period utilities.
                        # Then update expected values.
                
                        σ_1_t_plus_1 = σ_1[1, i+1]
                        κ_1_t_plus_1 = κ_1[:, i+1]
                
                        L_cohort_t_plus_1 = agg_labor_cohort(L_t_plus_1, σ_1_t_plus_1, κ_1_t_plus_1, N, R, C)
                
                        σ_0_t_plus_1 = σ_0[1, i+1]
                        κ_0_t_plus_1 = κ_0[:, i+1]
                
                
                        L_place_t_plus_1 = agg_labor_place(L_cohort_t_plus_1, σ_0_t_plus_1, κ_0_t_plus_1, N, R, C)
                
                        A_t_plus_1 = A[:, i+1]
                
                        w_t_plus_1 = wage(σ_0_t_plus_1, σ_1_t_plus_1, A_t_plus_1, κ_0_t_plus_1, κ_1_t_plus_1, N, R, C, L_place_t_plus_1, L_cohort_t_plus_1, L_t_plus_1)
                
                        r_bar_t_plus_1 = r_bar[:, i+1]
                        γ_t_plus_1 = γ[1, i+1]
                        η_t_plus_1 = η[1, i+1]
                
                
                
                        r_t_plus_1 = rent(w_t_plus_1, L_t_plus_1, r_bar_t_plus_1, η_t_plus_1, γ_t_plus_1, N, R, C)
                
                        B_t_plus_1 = B[:, i+1]
                
                        u_t_plus_1 = utility(w_t_plus_1, r_t_plus_1, B_t_plus_1, R, N, C, γ_t_plus_1)
                
                        s_t_plus_1 = s[:, i+1]
                        τ_t_plus_1 = τ[:, i+1]
                        V_t_plus_2 = V[:, i+2]
                        ν_t_plus_1 = ν[1, i+1]

                        V_t_plus_1 = expected_value_transition(u_t_plus_1, s_t_plus_1, ν_t_plus_1, τ_t_plus_1, R, N, C, V_t_plus_2)
                
                        V_new[:, i+1] = V_t_plus_1
                        
                        L[:, i+1] = L_t_plus_1  
                        w[:, i+1] = w_t_plus_1
                        r[:, i+1] = r_t_plus_1
                        
                        end
                end

                dif = maximum(abs.((V_new - V)./V))

                count = count + 1

                V = λ_2 * V_new + (1 - λ_2) * V


        end

        # γ[:, 1:(T-1)] is a 1×(T-1) matrix, γ[1, 1:(T-1)] is a (T-1) vector.
        real_wage = w ./ (repeat(kron(r, ones(C-1, 1)), R, 1) .^ (repeat(γ[:, 1:(T-1)], R*N*(C-1), 1)))



        # a bad feature of Julia is that it doesn't have lists in the language of R.
        # So, the output for this function should be matrix.
        # I add zeros as a "filler," because real wages and expected values/populations have different dimensions.

        real_wage = [real_wage; zeros(R*N, T-1)]

        output = [fill(count, N*R*C) fill(dif, N*R*C) V L real_wage]
        return output
end

path1 = transition_path(R, C, N, T, λ_2, L_in, V_in, V_in_2, s, ν, τ, M, σ_0, σ_1, κ_0, κ_1, A, r_bar, η, γ, B, tol, maxit)
V1 = path1[:, 3:2+T]
L1 = path1[:, (3+T):(2+2*T)]
real_wage1 = path1[1:N*R*(C-1), (3+2*T):(2+3*T-1)]

A_2 = zeros(N, T)
A_2[:, :] = A[:, :]
A_2[:, 4] = [8.0, 2.0]

path2 = transition_path(R, C, N, T, λ_2, L_in, V_in, V_in_2, s, ν, τ, M, σ_0, σ_1, κ_0, κ_1, A_2, r_bar, η, γ, B, tol, maxit)

V2 = path2[:, 3:3+(T-1)]
L2 = path2[:, 3+(T-1)+1:3+2*(T-1)]
real_wage2 = path2[1:N*R*(C-1), (3+2*T):(2+3*T-1)]
real_wage2

# compare expected values
# plot(0:10, V1[3, :1:11], label = "location 1 (baseline)")
# plot!(0:10, V1[C+3, :1:11], label = "location 2 (baseline)")
# plot!(0:10, V2[3, :1:11], label = "location 1 (counterfactual)")
# plot!(0:10, V2[C+3, :1:11], label = "location 2 (counterfactual)")
# xlabel!("period")
# ylabel!("expected value")
# savefig("../output/transition/prod_ex_val_age2.pdf")

# plot(0:10, V1[1, :1:11], label = "location 1 (baseline)")
# plot!(0:10, V1[C+1, :1:11], label = "location 2 (baseline)")
# plot!(0:10, V2[1, :1:11], label = "location 1 (counterfactual)")
# plot!(0:10, V2[C+1, :1:11], label = "location 2 (counterfactual)")
# xlabel!("period")
# ylabel!("expected value")
# savefig("../output/transition/prod_ex_val_age0.pdf")

value_cohort_path1_011 = zeros(C)
for i in 1:C
        value_cohort_path1_011[i] = V1[i, i]
end

value_cohort_path2_011 = zeros(C)
for i in 1:C
        value_cohort_path2_011[i] = V2[i, i]
end

value_cohort_path1_012 = zeros(C)
for i in 1:C
        value_cohort_path1_012[i] = V1[C+i, i]
end

value_cohort_path2_012 = zeros(C)
for i in 1:C
        value_cohort_path2_012[i] = V2[C+i, i]
end

# plot(0:(C-1), value_cohort_path1_011, label = "location 1 (baseline)")
# plot!(0:(C-1), value_cohort_path2_011, label = "location 1 (counterfactual)")
# plot!(0:(C-1), value_cohort_path1_012, label = "location 2 (baseline)")
# plot!(0:(C-1), value_cohort_path2_012, label = "location 2 (counterfactual)")
# xlabel!("age")
# ylabel!("expected value")
# savefig("../output/transition/expected_value_cohort_0.pdf")

# population
# plot(0:10, L1[3, :1:11], label = "location 1 (baseline)")
# plot!(0:10, L1[C+3, :1:11], label = "location 2 (baseline)")
# plot!(0:10, L2[3, :1:11], label = "location 1 (counterfactual)")
# plot!(0:10, L2[C+3, :1:11], label = "location 2 (counterfactual)")
# xlabel!("period")
# ylabel!("population")
# savefig("../output/transition/prod_pop_age2.pdf")

# plot(0:10, L1[3, :1:11], label = "age 2 in location 1 (baseline)")
# plot!(0:10, L1[C, :1:11], label = "age 7 in location 1 (baseline)")
# plot!(0:10, L2[3, :1:11], label = "age 2 in location 1 (counterfactual)")
# plot!(0:10, L2[C, :1:11], label = "age 7 in location 1 (counterfactual)")
# xlabel!("period")
# ylabel!("population")
# savefig("../output/transition/prod_pop_age2_7.pdf")

# real wage
# plot(0:10, real_wage1[2, :1:11], label = "location 1 (baseline)")
# plot!(0:10, real_wage1[(C-1)+2, :1:11], label = "location 2 (baseline)")
# plot!(0:10, real_wage2[2, :1:11], label = "location 1 (counterfactual)")
# plot!(0:10, real_wage2[(C-1)+2, :1:11], label = "location 2 (counterfactual)")
# xlabel!("period")
# ylabel!("real wage")
# savefig("../output/transition/prod_real_wage_age2.pdf")

# plot(0:10, real_wage1[1, :1:11], label = "location 1 (baseline)")
# plot!(0:10, real_wage1[(C-1)+1, :1:11], label = "location 2 (baseline)")
# plot!(0:10, real_wage2[1, :1:11], label = "location 1 (counterfactual)")
# plot!(0:10, real_wage2[(C-1)+1, :1:11], label = "location 2 (counterfactual)")
# xlabel!("period")
# ylabel!("real wage")

# plot(0:10, real_wage1[7, :1:11], label = "location 1 (baseline)")
# plot!(0:10, real_wage1[(C-1)+7, :1:11], label = "location 2 (baseline)")
# plot!(0:10, real_wage2[7, :1:11], label = "location 1 (counterfactual)")
# plot!(0:10, real_wage2[(C-1)+7, :1:11], label = "location 2 (counterfactual)")
# xlabel!("period")
# ylabel!("real wage")

# plot(0:10, real_wage1[2, :1:11], label = "location 1 (baseline)")
# plot!(0:10, real_wage2[2, :1:11], label = "location 1 (counterfactual)")
# plot!(0:10, real_wage1[7, :1:11], label = "location 1 (baseline)")
# plot!(0:10, real_wage2[7, :1:11], label = "location 1 (counterfactual)")


# plot(1:7, [real_wage1[1, 1], real_wage1[2, 2], real_wage1[3, 3], real_wage1[4, 4], real_wage1[5, 5], real_wage1[6, 6], real_wage1[7, 7]] )
# plot!(1:7, [real_wage2[1, 1], real_wage2[2, 2], real_wage2[3, 3], real_wage2[4, 4], real_wage2[5, 5], real_wage2[6, 6], real_wage2[7, 7]] )

# 011 means the cohort-race-place tuple
real_wage_cohort_path1_011 = zeros(C-1)
for i in 1:C-1
        real_wage_cohort_path1_011[i] = real_wage1[i, i+1]
end

real_wage_cohort_path2_011 = zeros(C-1)
for i in 1:C-1
        real_wage_cohort_path2_011[i] = real_wage2[i, i+1]
end

real_wage_cohort_path1_012 = zeros(C-1)
for i in 1:C-1
        real_wage_cohort_path1_012[i] = real_wage1[(C-1)+i, i+1]
end

real_wage_cohort_path2_012 = zeros(C-1)
for i in 1:C-1
        real_wage_cohort_path2_012[i] = real_wage2[(C-1)+i, i+1]
end

# plot(1:(C-1), real_wage_cohort_path1_011, label = "location 1 (baseline)")
# plot!(1:(C-1), real_wage_cohort_path2_011, label = "location 1 (counterfactual)")
# plot!(1:(C-1), real_wage_cohort_path1_012, label = "location 2 (baseline)")
# plot!(1:(C-1), real_wage_cohort_path2_012, label = "location 2 (counterfactual)")
# xlabel!("age")
# ylabel!("real wage")
# savefig("../output/transition/prod_real_wage_cohort_0.pdf")


# compute a steady state and a transition path for 50 places 

N = 50

# amenities 
B_period = fill(1.0, (R*N*C))

B = repeat(B_period, 1, T)

# survival probabilities
s_period = repeat([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0], R)
# the survival probability for the youngest race 1, ..., the survival probability for the second oldest race 1, ...
# the survival probability for the youngest race 2, ..., the survival probability for the second oldest race 2, ...
s = repeat(s_period, 1, T)


# migration costs 
# note that the oldest can no longer migrate

temp = fill(1.0, N, N)
temp = temp - Diagonal(fill(1.0, N))
τ_period = repeat(temp, R, C-1)
τ_period = reshape(τ_period', (N*R*N*(C-1)))

τ = repeat(τ_period, 1, T)

# place-specific shifter of rent 
r_bar_period = fill(0.1, N)

r_bar = repeat(r_bar_period, 1, T)

# place-specific elasticity of rent 
η_period = fill(0.5, N)

η = repeat(η_period, 1, T)

# productivity
A_period = Vector(51.0:100.0)

# the productivity in place 1, the productivity in place 2

A = repeat(A_period, 1, T)

# cohort-specific productivity 
κ_0_period = fill(1/7, N*(C-1))

κ_0 = repeat(κ_0_period, 1, T)

# race-specific productivity (within cohorts)

κ_1_period = fill(1/2, N*R*(C-1))
# the relative productivity of race1 within the oldest in place 1, ..., the relative productivity of race1 within the second youngest in place 1,
# the relative productivity of race1 within the oldest in place 2, ..., the relative productivity of race1 within the second youngest in place 2,
# the relative productivity of race2 within the oldest in place 1, ..., the relative productivity of race2 within the second youngest in place 1,
# the relative productivity of race2 within the oldest in place 2, ..., the relative productivity of race2 within the second youngest in place 2

κ_1 = repeat(κ_1_period, 1, T)

# fertility per cohort-race-place 
# I guess I need to restrict the value of α to get a steady state, 
α_period = repeat([0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0], R*N)
# the fertility of the second youngest race 1 in place 1, ..., the fertility of the oldest race 1 in place 1,
# the fertility of the second youngest race 2 in place 1, ..., the fertility of the oldest race 2 in place 1,
# the fertility of the second youngest race 2 in place 1, ..., the fertility of the oldest race 2 in place 1,
# the fertility of the second youngest race 2 in place 2, ..., the fertility of the oldest race 2 in place 2

α = repeat(α_period, 1, T)

# immgrants from abroad 
M_period = repeat([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], R*N)
# the immigrants of the youngest race 1 in place 1, ..., the immigrants of the oldest race 1 in place 1,
# the immigrants of the youngest race 1 in place 2, ..., the immigrants of the oldest race 1 in place 2,
# ...,
# the immigrants of the youngest race 2 in place 2, ..., the immigrants of the oldest race 2 in place 2.

M = repeat(M_period, 1, T)


# prepare for comparative transition paths

L_in = fill(10, (R*N*C, T))

LV_100 = steady_state(tol, maxit, N, C, R, ν[1, T], σ_0[1, T], σ_1[1, T], γ[1, T], B[:, T], s[:, T], τ[:, T], r_bar[:, T], η[:, T], A[:, T], κ_0[:, T], κ_1[:, T], α[:, T], L_in[:, T], λ)

V_mat_100_2 = reshape(LV_100[:, 2], C, R*N)'
L_mat_100_2 = reshape(LV_100[:, 1], C, R*N)'

L_mat_100 = reshape(LV_100[:, 1], C, R*N)'
V_mat_100 = reshape(LV_100[:, 2], C, R*N)'
real_wage_mat_100 = reshape(LV_100[1:N*R*(C-1), 3], C-1, R*N)'

maximum(abs.(L_mat_100[1, :] - L_mat_100[3, :]))

# I compute populations forward, given expected values.
# Then given populations, I update expected values.

λ_2 = 1.0

L_in = fill(10.0, (R*N*C), T)

V_in = fill(1.0, (R*N*C), T)
V_in_2 = fill(2.0, (R*N*C), T)

# the last period is the steady state
L_in[:, T] = LV_100[:, 1]
V_in[:, T] = LV_100[:, 2]
V_in_2[:, T] = LV_100[:, 2]

# the initial period (0) is the steady state
L_in[:, 1] = LV_100[:, 1]
V_in[:, 1] = LV_100[:, 2]
V_in_2[:, 1] = LV_100[:, 2]
# the "perceived" expected value before an unanticipated (MIT) shock happens in period 1

path1 = transition_path(R, C, N, T, λ_2, L_in, V_in, V_in_2, s, ν, τ, M, σ_0, σ_1, κ_0, κ_1, A, r_bar, η, γ, B, tol, maxit)
V1 = path1[:, 3:2+T]
L1 = path1[:, (3+T):(2+2*T)]
real_wage1 = path1[1:N*R*(C-1), (3+2*T):(2+3*T-1)]

A_2 = copy(A)
A_2[1, 4] = 150
A_2

path2 = transition_path(R, C, N, T, λ_2, L_in, V_in, V_in_2, s, ν, τ, M, σ_0, σ_1, κ_0, κ_1, A_2, r_bar, η, γ, B, tol, maxit)

V2 = path2[:, 3:3+(T-1)]
L2 = path2[:, 3+(T-1)+1:3+2*(T-1)]
real_wage2 = path2[1:N*R*(C-1), (3+2*T):(2+3*T-1)]
real_wage2