# Motoaki Takahashi 

# Julia 1.7.2

# March, April, August, September 2022

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

# suppose the economy converges to a steady state in 100 periods from 2010
# from period 1 (1940) to period 8 (2010), I feed fundamentals backed out from data 
# from period 9, I fix fundamentals and the economy converges to the steady state 
# in period 108 (8+100), the economy reaches to the steady state

T = 118

# the number of places 
N = 38

# years in the sample 
years = 1940:10:2010

# the number of cohorts in a given period
C = length(years)

# the number of races 
R = 2

# the number of years in data (1940, 1950, ..., 2010)
K = 8

# for parameters, columns generically represent periods,
# rows are races/places 

# migration elesticity
ν = fill(1.328, 1, T)

# elasticity of substitution between ages 
σ_0 = fill(1/0.4488, 1, T) 

# elasticity of substitution between races (within ages)
σ_1 = fill(1/0.118985, 1, T)

# Cobb-Douglas share of housing
γ = fill(0.25, 1, T)

# amenities 
amenity_data = CSV.read("../mig_wage/output/csv/amenity_matrix.csv")
# amenity_data = CSV.read("../mig_wage/output/csv/amenity_matrix_oridess.csv")

amenity_mat = convert(Matrix, amenity_data)
amenity_1940 = amenity_mat[:, 1] # suppose that amenities in 1940 are the same as those in 1950
amenity_mat2 = repeat(amenity_mat[:, 6], 1, (T-7)) # suppose that amenities since 2010 are those in 2000
B = hcat(amenity_1940, amenity_mat, amenity_mat2)

# columns are periods
# each row is structured as follows:
# the amenity for the youngest race 1 in place 1, ..., the amenity for the oldest race 1 in place 1, ...
# the amenity for the youngest race 1 in place 2, ..., the amenity for the oldest race 1 in place 2, ...
# the amenity for the youngest race 2 in place 1, ..., the amenity for the oldest race 2 in place 1, ...
# the amenity for the youngest race 2 in place 2, ..., the amenity for the oldest race 2 in place 2, ...



# survival probabilities
s_data = CSV.read("../life_table/output/s_matrix_q.csv")
s_mat = convert(Matrix, s_data)

s_mat2 = repeat(s_mat[:, 8], 1, (T-K)) # suppose that survival probabilities after 2010 are those in 2010

s = hcat(s_mat, s_mat2)
# columns are periods 
# each row is structured as:
# the survival probability for the youngest race 1, ..., the survival probability for the second oldest race 1, ...
# the survival probability for the youngest race 2, ..., the survival probability for the second oldest race 2, ...

# migration costs 
# note that the oldest can no longer migrate
τ_data = CSV.read("../mig_wage/output/csv/tau_matrix.csv")
# τ_data = CSV.read("../mig_wage/output/csv/tau_matrix_oridess.csv")

K_τ = 8
C_τ = 7

function make_τ_mat(τ_data, K_τ, C_τ)
        τ_mat = convert(Matrix, τ_data)
        τ_1 = zeros(R*C_τ*N*N, K_τ)

        for i in 1:K_τ
                temp = τ_mat[((i-1)R*N+1):(i*R*N) , :]'
                temp = reshape(temp, R*C_τ*N*N)
                τ_1[:, i] = temp
        end

        τ_2 = repeat(τ_1[:, K_τ], 1, (T-K_τ))

        τ = hcat(τ_1, τ_2)

        return τ
end

τ = make_τ_mat(τ_data, K_τ, C_τ)

# in τ_data, N×N migration cost matrices are filled as follows.
# columns: ages
# rows: (i) years in the outer tier, (ii) races in the inner tier

# in the matrix this code uses, columns are periods, and rows are structured as follows:


#  τ_period = [     0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
#                  1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0,
#                  0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
#                  1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0]

# the 1st row:
# the youngest race 1's migration cost from place 1 to place 1, the youngest race 1's migration cost from place 1 to place 2, ...,
# the second oldest race 1's migration cost from place 1 to place 1, the second oldest race 1's migration cost from place 1 to place 2, ...,
# ...,
# the last row:
# the youngest race 2's migration cost from place 2 to place 1, the youngest race 2's migration cost from place 2 to place 2, ...,
# the second oldest race 2's migration cost from place 2 to place 1, the second oldest race 2's migration cost from place 2 to place 2

# the outer tier: races
# the middle tier: ages 
# the lower tier: location×location

# place-specific shifter of rent 
r_bar_data = CSV.read("../rent_elas/output/csv/r_bar_2.csv")

r_bar_period = r_bar_data.r_bar_2

r_bar = repeat(r_bar_period, 1, T)

# place-specific elasticity of rent 
η_period = fill(0.4092, N)

η = repeat(η_period, 1, T)


# productivity
A_data = CSV.read("../elas_sub/output/csv/loc_prod_matrix.csv")

A_mat_1 =  convert(Matrix, A_data)

A_mat_2 = repeat(A_mat_1[:, K], 1, (T-K))


A = hcat(A_mat_1, A_mat_2)

A[:, 2] - A[:, 1]

# cohort-specific productivity 
κ_0_data = CSV.read("../elas_sub/output/csv/age_prod_matrix.csv")
κ_0_mat_1 = convert(Matrix, κ_0_data)

κ_0_mat_2 = repeat(κ_0_mat_1[:, K], 1, (T-K))


κ_0 = hcat(κ_0_mat_1, κ_0_mat_2)

# κ_0_period = [1/7, 1/7, 1/7, 1/7, 1/7, 1/7, 1/7,
#               1/7, 1/7, 1/7, 1/7, 1/7, 1/7, 1/7]
# the relative productivity of the second youngest in place 1, ..., the relative productivity of the oldest in place 1,
# the relative productivity of the second youngest in place 2, ..., the relative productivity of the oldest in place 2



# race-specific productivity (within cohorts)
κ_1_data = CSV.read("../elas_sub/output/csv/race_prod_matrix.csv")
κ_1_mat_1 = convert(Matrix, κ_1_data)

κ_1_mat_2 = repeat(κ_1_mat_1[:, K], 1, (T-K))

κ_1 = hcat(κ_1_mat_1, κ_1_mat_2)

# κ_1_period = [ 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2,
#                1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2,
#                1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2,
#                1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2]
# the relative productivity of race1 within the oldest in place 1, ..., the relative productivity of race1 within the second youngest in place 1,
# the relative productivity of race1 within the oldest in place 2, ..., the relative productivity of race1 within the second youngest in place 2,
# the relative productivity of race2 within the oldest in place 1, ..., the relative productivity of race2 within the second youngest in place 1,
# the relative productivity of race2 within the oldest in place 2, ..., the relative productivity of race2 within the second youngest in place 2



# fertility per cohort-race-place 
α_data = CSV.read("../fertility/output/csv/alpha_matrix.csv")
α_mat_1 = convert(Matrix, α_data)

# suppose that from period 9, age 20s (30) have the babies so that populations will be constant at the steady state 
# number of babies for nonblacks 
α_nb = 1/(s_mat2[1, 100] * s_mat2[2, 100])
α_b = 1/(s_mat2[8, 100] * s_mat2[9, 100])

# columns are periods.
# in rows:
# (i) races are in the upper tier,
# (ii) locations are in the middle tier,
# (iii) ages in the lower tier.

α_period_3 = vcat(repeat([0, 0, α_nb, 0, 0, 0, 0, 0], N), repeat([0, 0, α_b, 0, 0, 0, 0, 0], N))

α_mat_3 = repeat(α_period_3, 1, 100)

# 10 years of transition follow year 2010
# fertility in the transition period is between the fertility of 2010 and one of the steady state

α_mat_2 = zeros(R*N*C, 10)

for i in 1:10
        α_mat_2[:, i] = (i/10) * α_period_3 + ((10-i)/10) * α_mat_1[:,K]
end
α_mat_2


α = hcat(α_mat_1, α_mat_2, α_mat_3)

# I guess I need to restrict the value of α to get a steady state, 
# α_period = [    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
#                 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
#                 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
#                 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,]
# the fertility of the second youngest race 1 in place 1, ..., the fertility of the oldest race 1 in place 1,
# the fertility of the second youngest race 1 in place 2, ..., the fertility of the oldest race 1 in place 2, ???
# the fertility of the second youngest race 2 in place 1, ..., the fertility of the oldest race 2 in place 1,
# the fertility of the second youngest race 2 in place 2, ..., the fertility of the oldest race 2 in place 2



# immgrants from abroad 
# M_data = CSV.read("../migration_flow/state/output/csv/newcomer__pos_matrix.csv")
# M_mat_1 = convert(Matrix, M_data)

# I assume that since period 9 (2020), noone comes to the US from abroad, and noone leaves the US to abroad

# M_mat_2 = repeat(zeros(R*C*N) , 1, 100)

# M = hcat(M_mat_1, M_mat_2)

M = repeat(zeros(R*C*N) , 1, T)

# M_period = [    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
#                 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
#                 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
#                 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
# the immigrants of the youngest race 1 in place 1, ..., the immigrants of the oldest race 1 in place 1,
# the immigrants of the youngest race 1 in place 2, ..., the immigrants of the oldest race 1 in place 2,
# ...,
# the immigrants of the youngest race 2 in place 2, ..., the immigrants of the oldest race 2 in place 2.


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
#    temp = sum(w_mat .* L_mat[:, 1:(C-1)], dims = 2)
        temp = sum(w_mat .* L_mat[:, 2:C], dims = 2)
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

L_data = CSV.read("../migration_flow/state/output/csv/pop_q.csv")
L_mat_1 = convert(Matrix, L_data)
L_mat_2 = repeat(L_mat_1[:, K], 1, (T-K))

L_in = hcat(L_mat_1, L_mat_2)


# L_in = fill(10, (R*N*C, T))
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

# compute the steady state 

ss = steady_state(tol, maxit, N, C, R, ν[1, T], σ_0[1, T], σ_1[1, T], γ[1, T], B[:, T], s[:, T], τ[:, T], r_bar[:, T], η[:, T], A[:, T], κ_0[:, T], κ_1[:, T], α[:, T], L_in[:, K], λ)

(ss[:, 1] - L_in[:, T]) ./ L_in[:, T]
sum(L_in[:, T])
sum(L_in[:, T])

(sum(ss[:, 1]) - sum(L_in[:, T])) / sum(L_in[:, T])


V_mat_100 = reshape(ss[:, 2], C, R*N)'
L_mat_100 = reshape(ss[:, 1], C, R*N)'

real_wage_mat_100 = reshape(ss[1:N*R*(C-1), 3], C-1, R*N)'


# I compute populations forward, given expected values.
# Then given populations, I update expected values.

λ_2 = 1.0

# L_in = fill(10.0, (R*N*C), T)

V_in = fill(1.0, (R*N*C), T)
V_in_2 = fill(2.0, (R*N*C), T)

# the last period is the steady state
L_in[:, T] = ss[:, 1]
V_in[:, T] = ss[:, 2]
V_in_2[:, T] = ss[:, 2]

# the initial period (0) is year 1940 from data
# the "perceived" expected value before an unanticipated (MIT) shock happens in period 1


function transition_path(R, C, N, T, λ_2, L_in, V_in, V_in_2, s, ν, τ, M, σ_0, σ_1, κ_0, κ_1, A, r_bar, η, γ, B, tol, maxit, α)
        L = copy(L_in)
        V = copy(V_in) 
        V_new = copy(V_in_2)
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

        wage_mat = [w; zeros(R*N, T-1)]

        output = [fill(count, N*R*C) fill(dif, N*R*C) V L real_wage wage_mat]
        return output
end

path1 = transition_path(R, C, N, T, λ_2, L_in, V_in, V_in_2, s, ν, τ, M, σ_0, σ_1, κ_0, κ_1, A, r_bar, η, γ, B, tol, maxit)
V1 = path1[:, 3:2+T]
L1 = path1[:, (3+T):(2+2*T)]
real_wage1 = path1[1:N*R*(C-1), (3+2*T):(2+3*T-1)]
nom_wage1 =  path1[1:N*R*(C-1), (2+3*T):(1+4*T-1)]



L_in

V1

L1

sum(L1[1:N*C, :], dims = 1)
sum(L1[N*C+1:2*N*C, :], dims = 1) ./ (sum(L1[1:N*C, :], dims = 1) + sum(L1[N*C+1:2*N*C, :], dims = 1))

# the following function maps a transition path of nominal wages and populations to 
# nation-wide GDPs and per capita GDPs
function GDP(path_res)
        w_small = path_res[1:N*R*(C-1), (2+3*T):(1+4*T-1)]
        L = path_res[:, (3+T):(2+2*T)]

        w_big = zeros(R*N*C , (T-1))

        for i in 1:(T-1)
        for j in 1:R 
                for n in 1:N 
                        w_big[((j-1)*N*C+(n-1)*C+1):((j-1)*N*C+n*C) , i] = [0; w_small[((j-1)*N*(C-1)+(n-1)*(C-1)+1):((j-1)*N*(C-1)+n*(C-1)), i]]
                end
        end
        end

        income = L[:, 1:(T-1)] .* w_big

        US_income = sum(income, dims = 1)

        pop = sum(L[:, 1:(T-1)], dims = 1)

        US_per_capita_income = US_income ./ pop

        output = [US_income; US_per_capita_income]
        return output
end

function local_GDP(path_res, N, T)
        w_small = path_res[1:N*R*(C-1), (2+3*T):(1+4*T-1)]
        L = path_res[:, (3+T):(2+2*T)]

        w_big = zeros(R*N*C , (T-1))

        for i in 1:(T-1)
        for j in 1:R 
                for n in 1:N 
                        w_big[((j-1)*N*C+(n-1)*C+1):((j-1)*N*C+n*C) , i] = [0; w_small[((j-1)*N*(C-1)+(n-1)*(C-1)+1):((j-1)*N*(C-1)+n*(C-1)), i]]
                end
        end
        end

        income = L[:, 1:(T-1)] .* w_big

        # compute location-level income 
        
        local_gdp = zeros(N, T-1)

        for i in 1:R  
                for n in 1:N 
                        local_gdp[n, :] = sum(income[((i-1)*N*C+(n-1)*C+1):((i-1)*N*C+n*C) , :], dims = 1)
                end
        end

        return local_gdp

end

gdp1 = GDP(path1)
local_gdp1 = local_GDP(path1, N, T)

# the following function maps a transition path of populations to the shares of African Americans and
# the others in the South

# read the concordance table mapping geographic units to the South 

south_dt = CSV.read("../geo_unit/output/my_geo_south.csv")

north_or_south = south_dt.south
# takes 0 if the location is in the North,
# takes 1 if the location is in the South.

south_dt = [south_dt 1:N]

south_ind_df = filter(row -> row.south == 1, south_dt)
south_ind = south_ind_df.x1

north_ind_df = filter(row -> row.south == 0, south_dt)
north_ind = north_ind_df.x1

N_S = length(south_ind)
N_N = length(north_ind)

function share_south(path_res, north_or_south)
        L = path_res[:, (3+T):(2+2*T)]

        south_pop = zeros(2*R, T)
        # columns are periods 
        # rows are:
        # (i) races are in the outer tier,
        # (ii) north/south are in the inner tier

        for t in 1:T 
                for i in 1:R 
                        for n in 1:N 
                                northorsouth = north_or_south[n]
                                temp = L[((i-1)*N*C+(n-1)*C+1):((i-1)*N*C+n*C), t]
                                south_pop[(i-1)*2+1+northorsouth , t] = south_pop[(i-1)*2+1+northorsouth , t] + sum(temp)
                        end
                end
        end

        south_share_b = south_pop[4, :] ./ (south_pop[3, :] + south_pop[4, :])
        south_share_nb = south_pop[2, :] ./ (south_pop[1, :] + south_pop[2, :])
        output = [south_share_nb'; south_share_b']

        return output
end

share_south1 = share_south(path1, north_or_south)


# no migration between the North and the South from 1940 to 1960
τ_nonorthsouth4060_data = CSV.read("../mig_wage/output/csv/tau_nonorthsouth4060_matrix.csv")

τ_nonorthsouth4060 = make_τ_mat(τ_nonorthsouth4060_data, K_τ, C_τ)

ss_nonorthsouth4060 = steady_state(tol, maxit, N, C, R, ν[1, T], σ_0[1, T], σ_1[1, T], γ[1, T], B[:, T], s[:, T], τ_nonorthsouth4060[:, T], r_bar[:, T], η[:, T], A[:, T], κ_0[:, T], κ_1[:, T], α[:, T], L_in[:, K], λ)

V_nonorthsouth4060_mat_100 = reshape(ss_nonorthsouth4060[:, 2], C, R*N)'
L_nonorthsouth4060_mat_100 = reshape(ss_nonorthsouth4060[:, 1], C, R*N)'

V_in = fill(1.0, (R*N*C), T)
V_in_2 = fill(2.0, (R*N*C), T)

# the last period is the steady state
L_in[:, T] = ss_nonorthsouth4060[:, 1]
V_in[:, T] = ss_nonorthsouth4060[:, 2]
V_in_2[:, T] = ss_nonorthsouth4060[:, 2]

path_nonorthsouth4060 = transition_path(R, C, N, T, λ_2, L_in, V_in, V_in_2, s, ν, τ_nonorthsouth4060, M, σ_0, σ_1, κ_0, κ_1, A, r_bar, η, γ, B, tol, maxit)
V_nonorthsouth4060 = path_nonorthsouth4060[:, 3:2+T]
L_nonorthsouth4060 = path_nonorthsouth4060[:, (3+T):(2+2*T)]
real_wage_nonorthsouth4060 = path_nonorthsouth4060[1:N*R*(C-1), (3+2*T):(2+3*T-1)]
nom_wage_nonorthsouth4060 =  path_nonorthsouth4060[1:N*R*(C-1), (2+3*T):(1+4*T-1)]


gdp_nonorthsouth4060 = GDP(path_nonorthsouth4060)

gdp_nonorthsouth4060 ./ gdp1

V_nonorthsouth4060 ./ V1

share_south_nonorthsouth4060 = share_south(path_nonorthsouth4060, north_or_south)

# no migration between the North and the South from 1940 to 1960 only for African Americans 
τ_nonorthsouth4060_b_data = CSV.read("../mig_wage/output/csv/tau_nonorthsouth4060_b_matrix.csv")

τ_nonorthsouth4060_b = make_τ_mat(τ_nonorthsouth4060_b_data, K_τ, C_τ)

ss_nonorthsouth4060_b = steady_state(tol, maxit, N, C, R, ν[1, T], σ_0[1, T], σ_1[1, T], γ[1, T], B[:, T], s[:, T], τ_nonorthsouth4060_b[:, T], r_bar[:, T], η[:, T], A[:, T], κ_0[:, T], κ_1[:, T], α[:, T], L_in[:, K], λ)

V_nonorthsouth4060_b_mat_100 = reshape(ss_nonorthsouth4060_b[:, 2], C, R*N)'
L_nonorthsouth4060_b_mat_100 = reshape(ss_nonorthsouth4060_b[:, 1], C, R*N)'

V_in = fill(1.0, (R*N*C), T)
V_in_2 = fill(2.0, (R*N*C), T)

# the last period is the steady state
L_in[:, T] = ss_nonorthsouth4060_b[:, 1]
V_in[:, T] = ss_nonorthsouth4060_b[:, 2]
V_in_2[:, T] = ss_nonorthsouth4060_b[:, 2]

path_nonorthsouth4060_b = transition_path(R, C, N, T, λ_2, L_in, V_in, V_in_2, s, ν, τ_nonorthsouth4060_b, M, σ_0, σ_1, κ_0, κ_1, A, r_bar, η, γ, B, tol, maxit)
V_nonorthsouth4060_b = path_nonorthsouth4060_b[:, 3:2+T]
L_nonorthsouth4060_b = path_nonorthsouth4060_b[:, (3+T):(2+2*T)]
real_wage_nonorthsouth4060_b = path_nonorthsouth4060_b[1:N*R*(C-1), (3+2*T):(2+3*T-1)]
nom_wage_nonorthsouth4060_b =  path_nonorthsouth4060_b[1:N*R*(C-1), (2+3*T):(1+4*T-1)]

gdp_nonorthsouth4060_b = GDP(path_nonorthsouth4060_b)

gdp_ratio_nonorthsouth4060_b = gdp_nonorthsouth4060_b ./ gdp1

plot(years, gdp_ratio_nonorthsouth4060_b[1, 1:K], xlabel = "year",
        ylabel ="GDP relative to the baseline", legend = false)
savefig("output/figures/gdp_nonorthsouth4060_b.pdf")

local_gdp_nonorthsouth4060_b = local_GDP(path_nonorthsouth4060_b, N, T)

local_gdp_nonorthsouth4060_b_ratio = local_gdp_nonorthsouth4060_b ./ local_gdp1

plot(years, local_gdp_nonorthsouth4060_b_ratio[21, 1:K], xlabel = "year", ylabel = "GDP relative to the baseline", label = "Mississippi",
 c = :black)
plot!(years, local_gdp_nonorthsouth4060_b_ratio[11, 1:K], label = "Illinois", c = :black, line = :dash)
savefig("output/figures/gdp_nonorthsouth4060_b_MS_IL.pdf")


# collect values of the youngest
V_nonorthsouth4060_b

function get_V_youngest(V)
        output = zeros(R*N, K)
        for k in 1:K
                for i in 1:R
                        for j in 1:N
                                output[(i-1)*N+j, k] = V[(i-1)*N*C+(j-1)*C+1, k]
                        end
                end
        end
        
        return output

end


V_nonorthsouth4060_b_youngest = get_V_youngest(V_nonorthsouth4060_b)
V1_youngest = get_V_youngest(V1)

V_nonorthsouth4060_b_youngest_ratio = V_nonorthsouth4060_b_youngest ./ V1_youngest

V_nonorthsouth4060_b_youngest_ratio[N+1:R*N,:]

V_nonorthsouth4060_b_youngest_ratio[21,:]
V_nonorthsouth4060_b_youngest_ratio[11,:]

plot(years - fill(10.0, K), V_nonorthsouth4060_b_youngest_ratio[N+21, 1:K], xlabel = "cohort",
ylabel = "Expected value at the beggining of life", label = "Mississippi", legend = :bottomright, c=:black)
plot!(years - fill(10.0, K), V_nonorthsouth4060_b_youngest_ratio[N+11, 1:K], label = "Illinois", c=:black,
line=:dash)
savefig("output/figures/cohort_V_nonorthsouth4060_b_MS_IL.pdf")

plot(years - fill(10.0, K), V_nonorthsouth4060_b_youngest_ratio[21, 1:K], xlabel = "cohort",
ylabel = "Expected value at the beggining of life", label = "Mississippi", legend = :bottomright, c=:black)
plot!(years - fill(10.0, K), V_nonorthsouth4060_b_youngest_ratio[11, 1:K], label = "Illinois", c=:black,
line=:dash)
savefig("output/figures/cohort_V_nonorthsouth4060_b_MS_IL_for_o.pdf")

# compute the average value of the youngest for the South and the North in each year

function region_youngest_V(V_youngest, N, R, K, north_or_south, N_S, N_N)
        output = zeros(2*R, K)

        for i in 1:R
        for n in 1:N 
                if north_or_south[n] == 1 # if the locaiton is in the South
                        output[(i-1)*2+1, :] = output[(i-1)*2+1, :] + V_youngest[(i-1)*N+n, 1:K] ./ N_S
                else # if the location is in the North
                        output[(i-1)*2+2, :] = output[(i-1)*2+2, :] + V_youngest[(i-1)*N+n, 1:K] ./ N_N
                end
        end
        end
        return output
end

V_youngest_NS_nonorthsouth4060_b = region_youngest_V(V_nonorthsouth4060_b, N, R, K, north_or_south, N_S, N_N)
V_youngest_NS_1 = region_youngest_V(V1, N, R, K, north_or_south, N_S, N_N)

V_youngest_ratio_nonorthsouth4060_b = V_youngest_NS_nonorthsouth4060_b ./ V_youngest_NS_1

plot(years, V_youngest_ratio_nonorthsouth4060_b[1, :])
plot!(years, V_youngest_ratio_nonorthsouth4060_b[2, :])
plot!(years, V_youngest_ratio_nonorthsouth4060_b[3, :])
plot!(years, V_youngest_ratio_nonorthsouth4060_b[4, :])

# compare standard deviations in expected values within races 
plot(years - fill(10.0, K), std(V1_youngest[N+1:2*N, :], dims=1)', xlabel = "cohort",
ylabel = "standard deviation of expected values", label = "baseline")
plot!(years - fill(10.0, K), std(V_nonorthsouth4060_b_youngest[N+1:2*N, :], dims=1)',
        label = "no Black migration across the North and the South")

plot(years - fill(10.0, K), std(V1_youngest[1:N, :], dims=1)', xlabel = "cohort",
        ylabel = "standard deviation of expected values", label = "baseline")
plot!(years - fill(10.0, K), std(V_nonorthsouth4060_b_youngest[1:N, :], dims=1)',
        label = "no Black migration across the North and the South")






#  no migration between the North and the South from 1940 to 1960 only for nonblacks 
τ_nonorthsouth4060_nb_data = CSV.read("../mig_wage/output/csv/tau_nonorthsouth4060_nb_matrix.csv")

τ_nonorthsouth4060_nb = make_τ_mat(τ_nonorthsouth4060_nb_data, K_τ, C_τ)

ss_nonorthsouth4060_nb = steady_state(tol, maxit, N, C, R, ν[1, T], σ_0[1, T], σ_1[1, T], γ[1, T], B[:, T], s[:, T], τ_nonorthsouth4060_nb[:, T], r_bar[:, T], η[:, T], A[:, T], κ_0[:, T], κ_1[:, T], α[:, T], L_in[:, K], λ)

V_nonorthsouth4060_nb_mat_100 = reshape(ss_nonorthsouth4060_nb[:, 2], C, R*N)'
L_nonorthsouth4060_nb_mat_100 = reshape(ss_nonorthsouth4060_nb[:, 1], C, R*N)'

V_in = fill(1.0, (R*N*C), T)
V_in_2 = fill(2.0, (R*N*C), T)

# the last period is the steady state
L_in[:, T] = ss_nonorthsouth4060_nb[:, 1]
V_in[:, T] = ss_nonorthsouth4060_nb[:, 2]
V_in_2[:, T] = ss_nonorthsouth4060_nb[:, 2]

path_nonorthsouth4060_nb = transition_path(R, C, N, T, λ_2, L_in, V_in, V_in_2, s, ν, τ_nonorthsouth4060_nb, M, σ_0, σ_1, κ_0, κ_1, A, r_bar, η, γ, B, tol, maxit)
V_nonorthsouth4060_nb = path_nonorthsouth4060_nb[:, 3:2+T]
L_nonorthsouth4060_nb = path_nonorthsouth4060_nb[:, (3+T):(2+2*T)]
real_wage_nonorthsouth4060_nb = path_nonorthsouth4060_nb[1:N*R*(C-1), (3+2*T):(2+3*T-1)]
nom_wage_nonorthsouth4060_nb =  path_nonorthsouth4060_nb[1:N*R*(C-1), (2+3*T):(1+4*T-1)]

gdp_nonorthsouth4060_nb = GDP(path_nonorthsouth4060_nb)

gdp_ratio_nonorthsouth4060_nb = gdp_nonorthsouth4060_nb ./ gdp1

plot(years, gdp_ratio_nonorthsouth4060_b[1, 1:K], xlabel = "year",
        ylabel ="US GDP relative to the baseline", label = "No migration of African Americans", c=:black)
plot!(years, gdp_ratio_nonorthsouth4060_nb[1, 1:K], label = "No migration of the Others", c=:black, line = :dash)

savefig("output/figures/gdp_nonorthsouth4060_b_vs_nb.pdf")



V_nonorthsouth4060_nb_youngest = get_V_youngest(V_nonorthsouth4060_nb)

V_nonorthsouth4060_nb_youngest_ratio = V_nonorthsouth4060_nb_youngest ./ V1_youngest

plot(years - fill(10.0, K), V_nonorthsouth4060_nb_youngest_ratio[21, 1:K], xlabel = "cohort",
ylabel = "Expected value at the beggining of life", label = "Mississippi", legend = :bottomright, c=:black)
plot!(years - fill(10.0, K), V_nonorthsouth4060_nb_youngest_ratio[11, 1:K], label = "Illinois", c=:black, line = :dash)
savefig("output/figures/cohort_V_nonorthsouth4060_nb_MS_IL.pdf")

local_gdp_nonorthsouth4060_nb = local_GDP(path_nonorthsouth4060_nb, N, T)

local_gdp_nonorthsouth4060_nb_ratio = local_gdp_nonorthsouth4060_nb ./ local_gdp1

plot(years, local_gdp_nonorthsouth4060_nb_ratio[21, 1:K], xlabel = "year", ylabel = "GDP relative to the baseline", 
        label = "Mississippi", c = :black)
plot!(years, local_gdp_nonorthsouth4060_nb_ratio[11, 1:K], xlabel = "year", ylabel = "GDP relative to the baseline", 
        label = "Illinois", c = :black, line = :dash)
savefig("output/figures/gdp_nonorthsouth4060_nb_MS_IL.pdf")

# counterfactual economy in which African Americans have the others' migration costs (suggested by Jonathan)

# τ is the factual migration costs 

τ_nb_to_b = zeros(size(τ))
N*N*R*(K-1)

τ_nb_to_b[1:(C_τ*N*N), :] = copy(τ[1:(C_τ*N*N), :])
τ_nb_to_b[(C_τ*N*N+1):(R*C_τ*N*N), :] = copy(τ[1:(C_τ*N*N), :])

ss_τ_nb_to_b = steady_state(tol, maxit, N, C, R, ν[1, T], σ_0[1, T], σ_1[1, T], γ[1, T], B[:, T], s[:, T], τ_nb_to_b[:, T], r_bar[:, T], η[:, T], A[:, T], κ_0[:, T], κ_1[:, T], α[:, T], L_in[:, K], λ)

V_τ_nb_to_b = reshape(ss_τ_nb_to_b[:, 2], C, R*N)'
L_τ_nb_to_b = reshape(ss_τ_nb_to_b[:, 1], C, R*N)'

V_in = fill(1.0, (R*N*C), T)
V_in_2 = fill(2.0, (R*N*C), T)

# the last period is the steady state
L_in[:, T] = ss_τ_nb_to_b[:, 1]
V_in[:, T] = ss_τ_nb_to_b[:, 2]
V_in_2[:, T] = ss_τ_nb_to_b[:, 2]

path_τ_nb_to_b = transition_path(R, C, N, T, λ_2, L_in, V_in, V_in_2, s, ν, τ_nb_to_b, M, σ_0, σ_1, κ_0, κ_1, A, r_bar, η, γ, B, tol, maxit)
V_τ_nb_to_b = path_τ_nb_to_b[:, 3:2+T]
L_τ_nb_to_b = path_τ_nb_to_b[:, (3+T):(2+2*T)]
real_wage_τ_nb_to_b = path_τ_nb_to_b[1:N*R*(C-1), (3+2*T):(2+3*T-1)]
nom_wage_τ_nb_to_b =  path_τ_nb_to_b[1:N*R*(C-1), (2+3*T):(1+4*T-1)]

gdp_τ_nb_to_b = GDP(path_τ_nb_to_b)

gdp_ratio_τ_nb_to_b = gdp_τ_nb_to_b ./ gdp1

V_τ_nb_to_b_youngest = get_V_youngest(V_τ_nb_to_b)

V_τ_nb_to_b_youngest_ratio = V_τ_nb_to_b_youngest ./ V1_youngest

plot(years - fill(10.0, K), V_τ_nb_to_b_youngest_ratio[N+21, 1:K], xlabel = "cohort",
ylabel = "Expected value at the beggining of life", label = "Mississippi", legend = :bottomright, c=:black)
plot!(years - fill(10.0, K), V_τ_nb_to_b_youngest_ratio[N+11, 1:K], label = "Illinois", c=:black, line = :dash)
# savefig("output/figures/cohort_V_tau_nb_to_b_youngest_ratio.pdf")

plot(years, gdp_ratio_τ_nb_to_b[1, 1:K], xlabel = "year",
        ylabel ="US GDP relative to the baseline", label = "African Americans have the others' migration costs", c=:black)

# savefig("output/figures/gdp_tau_nb_to_b.pdf")

plot(years, gdp_ratio_τ_nb_to_b[2, 1:K], xlabel = "year",
        ylabel ="US per capita GDP relative to the baseline", label = "African Americans have the others' migration costs", c=:black)

# savefig("output/figures/per_capita_gdp_tau_nb_to_b.pdf")

# African Americans have the others migration costs from 1940 to 1960

τ_nb_to_b_1940_1960 = copy(τ)

# African Americans' migration costs are as the others' migration costs from 1940 to 1960
τ_nb_to_b_1940_1960[(C_τ*N*N+1):(R*C_τ*N*N), 1:3] = copy(τ[1:(C_τ*N*N), 1:3])

ss_τ_nb_to_b_1940_1960 = steady_state(tol, maxit, N, C, R, ν[1, T], σ_0[1, T], σ_1[1, T], γ[1, T], B[:, T], s[:, T], τ_nb_to_b_1940_1960[:, T], r_bar[:, T], η[:, T], A[:, T], κ_0[:, T], κ_1[:, T], α[:, T], L_in[:, K], λ)

V_τ_nb_to_b_1940_1960 = reshape(ss_τ_nb_to_b_1940_1960[:, 2], C, R*N)'
L_τ_nb_to_b_1940_1960 = reshape(ss_τ_nb_to_b_1940_1960[:, 1], C, R*N)'

V_in = fill(1.0, (R*N*C), T)
V_in_2 = fill(2.0, (R*N*C), T)

# the last period is the steady state
L_in[:, T] = ss_τ_nb_to_b_1940_1960[:, 1]
V_in[:, T] = ss_τ_nb_to_b_1940_1960[:, 2]
V_in_2[:, T] = ss_τ_nb_to_b_1940_1960[:, 2]

path_τ_nb_to_b_1940_1960 = transition_path(R, C, N, T, λ_2, L_in, V_in, V_in_2, s, ν, τ_nb_to_b_1940_1960, M, σ_0, σ_1, κ_0, κ_1, A, r_bar, η, γ, B, tol, maxit)
V_τ_nb_to_b_1940_1960 = path_τ_nb_to_b_1940_1960[:, 3:2+T]
L_τ_nb_to_b_1940_1960 = path_τ_nb_to_b_1940_1960[:, (3+T):(2+2*T)]
real_wage_τ_nb_to_b_1940_1960 = path_τ_nb_to_b_1940_1960[1:N*R*(C-1), (3+2*T):(2+3*T-1)]
nom_wage_τ_nb_to_b_1940_1960 =  path_τ_nb_to_b_1940_1960[1:N*R*(C-1), (2+3*T):(1+4*T-1)]

gdp_τ_nb_to_b_1940_1960 = GDP(path_τ_nb_to_b_1940_1960)

gdp_ratio_τ_nb_to_b_1940_1960 = gdp_τ_nb_to_b_1940_1960 ./ gdp1

V_τ_nb_to_b_1940_1960_youngest = get_V_youngest(V_τ_nb_to_b_1940_1960)

V_τ_nb_to_b_1940_1960_youngest_ratio = V_τ_nb_to_b_1940_1960_youngest ./ V1_youngest

plot(years, gdp_ratio_τ_nb_to_b_1940_1960[1, 1:K], xlabel = "year",
        ylabel ="US GDP relative to the baseline", label = "1940-1960", c=:black, legend = :bottomright)
plot!(years, gdp_ratio_τ_nb_to_b[1, 1:K], label = "forever", c=:black, line=:dash)


savefig("output/figures/gdp_nb_to_b_1940_1960_vs_forever.pdf")

plot(years, gdp_ratio_τ_nb_to_b_1940_1960[2, 1:K], xlabel = "year",
        ylabel ="US per capita GDP relative to the baseline", label = "1940-1960", c=:black, legend = :bottomright)
plot!(years, gdp_ratio_τ_nb_to_b[2, 1:K], label = "forever", c=:black, line=:dash)


savefig("output/figures/per_capita_gdp_nb_to_b_1940_1960_vs_forever.pdf")

plot(years - fill(10.0, K), V_τ_nb_to_b_youngest_ratio[N+21, 1:K], xlabel = "cohort",
ylabel = "Expected value at the beggining of life", label = "Mississippi: forever", legend = :right, c=:black, line = :solid)
plot!(years - fill(10.0, K), V_τ_nb_to_b_youngest_ratio[N+11, 1:K], label = "Illinois: forever", c=:black, line = :dash)

plot!(years - fill(10.0, K), V_τ_nb_to_b_1940_1960_youngest_ratio[N+21, 1:K], label = "Mississippi: 1940-1960", c=:green, line = :solid)
plot!(years - fill(10.0, K), V_τ_nb_to_b_1940_1960_youngest_ratio[N+11, 1:K], label = "Illinois: 1940-1960", c=:green, line = :dash)

savefig("output/figures/expected_value_nb_to_b_1940_1960_vs_forever.pdf")


# African Americans cannot move from the South to the North from 1940 to 1960
τ_no_s_to_n_b_data = CSV.read("../mig_wage/output/csv/tau_no_south_to_north_4060_b_matrix.csv")

τ_no_s_to_n_b = make_τ_mat(τ_no_s_to_n_b_data, K_τ, C_τ)

ss_no_s_to_n_b = steady_state(tol, maxit, N, C, R, ν[1, T], σ_0[1, T], σ_1[1, T], γ[1, T], B[:, T], s[:, T], τ_no_s_to_n_b[:, T], r_bar[:, T], η[:, T], A[:, T], κ_0[:, T], κ_1[:, T], α[:, T], L_in[:, K], λ)

V_no_s_to_n_b_mat_100 = reshape(ss_no_s_to_n_b[:, 2], C, R*N)'
L_no_s_to_n_b_mat_100 = reshape(ss_no_s_to_n_b[:, 1], C, R*N)'

V_in = fill(1.0, (R*N*C), T)
V_in_2 = fill(2.0, (R*N*C), T)

# the last period is the steady state
L_in[:, T] = ss_no_s_to_n_b[:, 1]
V_in[:, T] = ss_no_s_to_n_b[:, 2]
V_in_2[:, T] = ss_no_s_to_n_b[:, 2]

path_no_s_to_n_b = transition_path(R, C, N, T, λ_2, L_in, V_in, V_in_2, s, ν, τ_no_s_to_n_b, M, σ_0, σ_1, κ_0, κ_1, A, r_bar, η, γ, B, tol, maxit)
V_no_s_to_n_b = path_no_s_to_n_b[:, 3:2+T]
L_no_s_to_n_b = path_no_s_to_n_b[:, (3+T):(2+2*T)]
real_wage_no_s_to_n_b = path_no_s_to_n_b[1:N*R*(C-1), (3+2*T):(2+3*T-1)]
nom_wage_no_s_to_n_b =  path_no_s_to_n_b[1:N*R*(C-1), (2+3*T):(1+4*T-1)]

gdp_no_s_to_n_b = GDP(path_no_s_to_n_b)

gdp_ratio_no_s_to_n_b = gdp_no_s_to_n_b ./ gdp1

V_no_s_to_n_b_youngest = get_V_youngest(V_no_s_to_n_b)

V_no_s_to_n_b_youngest_ratio = V_no_s_to_n_b_youngest ./ V1_youngest

plot(years, gdp_ratio_no_s_to_n_b[1, 1:K], xlabel = "year",
        ylabel ="US GDP relative to the baseline", label = "No migraiton from the South to the North", c=:black, legend = :bottomright,
        ylims = [0.984, 1.005])
plot!(years, gdp_ratio_nonorthsouth4060_b[1, 1:K], xlabel = "year",
        ylabel ="US GDP relative to the baseline", label = "No migraiton across the North and the South", c=:black,
        line = :dash)

savefig("output/figures/gdp_no_mig_4060_b_north_south_bilateral_vs_unilateral.pdf")

plot(years - fill(10.0, K), V_no_s_to_n_b_youngest_ratio[N+21, 1:K], label = "No migraiton from the South to the North",
c=:black, line = :solid, legend = :bottomright)
plot!(years - fill(10.0, K), V_nonorthsouth4060_b_youngest_ratio[N+21, 1:K], label = "No migraiton across the North and the South", c=:black, line = :dash)

savefig("output/figures/expected_value_no_mig_4060_b_north_south_bilateral_vs_unilateral.pdf")

plot(years, gdp_ratio_no_s_to_n_b[2, 1:K], xlabel = "year",
        ylabel ="US per capita GDP relative to the baseline", label = "No migraiton from the South to the North",
        c=:black, legend = :topleft,
        ylims = [0.984, 1.005])
plot!(years, gdp_ratio_nonorthsouth4060_b[1, 1:K], xlabel = "year",
        ylabel ="US per capita GDP relative to the baseline", label = "No migraiton across the North and the South", c=:black,
        line = :dash)





# counterfactual economy in which migration costs are zero 
τ_path2 = zeros(R*C_τ*N*N, T)

ss_path2 = steady_state(tol, maxit, N, C, R, ν[1, T], σ_0[1, T], σ_1[1, T], γ[1, T], B[:, T], s[:, T], τ_path2[:, T], r_bar[:, T], η[:, T], A[:, T], κ_0[:, T], κ_1[:, T], α[:, T], L_in[:, K], λ)

# the last period is the steady state
L_in_path2 = copy(L_in)
L_in_path2[:, T] = ss_path2[:, 1]

V_in_path2 = copy(V_in)
V_in_path2[:, T] = ss_path2[:, 2]

V_in_path2_2 = copy(V_in_2)
V_in_path2_2[:, T] = ss_path2[:, 2]

path2 = transition_path(R, C, N, T, λ_2, L_in_path2, V_in_path2, V_in_path2_2, s, ν, τ_path2, M, σ_0, σ_1, κ_0, κ_1, A, r_bar, η, γ, B, tol, maxit)

V2 = path2[:, 3:2+T]
L2 = path2[:, (3+T):(2+2*T)]
real_wage2 = path2[1:N*R*(C-1), (3+2*T):(2+3*T-1)]
nom_wage2 =  path2[1:N*R*(C-1), (2+3*T):(1+4*T-1)]

gdp2 = GDP(path2)

gdp2 ./ gdp1

# migration autarky from 1940 to 1960
τ_autarky4060_data = CSV.read("../mig_wage/output/csv/tau_autarky4060_matrix.csv")

τ_autarky4060 = make_τ_mat(τ_autarky4060_data, K_τ, C_τ)

ss_autarky4060 = steady_state(tol, maxit, N, C, R, ν[1, T], σ_0[1, T], σ_1[1, T], γ[1, T], B[:, T], s[:, T], τ_autarky4060[:, T], r_bar[:, T], η[:, T], A[:, T], κ_0[:, T], κ_1[:, T], α[:, T], L_in[:, K], λ)

# the last period is the steady state
L_in_autarky4060 = copy(L_in)
L_in_autarky4060[:, T] = ss_autarky4060[:, 1]

V_in_autarky4060 = copy(V_in)
V_in_autarky4060[:, T] = ss_autarky4060[:, 2]

V_in_autarky4060_2 = copy(V_in_2)
V_in_autarky4060_2[:, T] = ss_autarky4060[:, 2]

path_autarky4060 = transition_path(R, C, N, T, λ_2, L_in_autarky4060, V_in_autarky4060, V_in_autarky4060_2, s, ν, τ_autarky4060, M, σ_0, σ_1, κ_0, κ_1, A, r_bar, η, γ, B, tol, maxit)

V_autarky4060 = path_autarky4060[:, 3:2+T]
L_autarky4060 = path_autarky4060[:, (3+T):(2+2*T)]
real_wage_autarky4060 = path_autarky4060[1:N*R*(C-1), (3+2*T):(2+3*T-1)]
nom_wage_autarky4060 =  path_autarky4060[1:N*R*(C-1), (2+3*T):(1+4*T-1)]

gdp_autarky4060 = GDP(path_autarky4060)

gdp_autarky4060 ./ gdp1

# insert migration costs as of 2010 to every year 
τ_always2010 = repeat(τ[:,8], 1, T)

ss_always2010 = steady_state(tol, maxit, N, C, R, ν[1, T], σ_0[1, T], σ_1[1, T], γ[1, T], B[:, T], s[:, T], τ_always2010[:, T], r_bar[:, T], η[:, T], A[:, T], κ_0[:, T], κ_1[:, T], α[:, T], L_in[:, K], λ)

# the last period is the steady state
L_in_always2010 = copy(L_in)
L_in_always2010[:, T] = ss_always2010[:, 1]

V_in_always2010 = copy(V_in)
V_in_always2010[:, T] = ss_always2010[:, 2]

V_in_always2010_2 = copy(V_in_2)
V_in_always2010_2[:, T] = ss_always2010[:, 2]

path_always2010 = transition_path(R, C, N, T, λ_2, L_in_always2010, V_in_always2010, V_in_always2010_2, s, ν, τ_always2010, M, σ_0, σ_1, κ_0, κ_1, A, r_bar, η, γ, B, tol, maxit)

V_always2010 = path_always2010[:, 3:2+T]
L_always2010 = path_always2010[:, (3+T):(2+2*T)]
real_wage_always2010 = path_always2010[1:N*R*(C-1), (3+2*T):(2+3*T-1)]
nom_wage_always2010 =  path_always2010[1:N*R*(C-1), (2+3*T):(1+4*T-1)]

gdp_always2010 = GDP(path_always2010)

gdp_always2010 ./ gdp1

# African Americans have the others' migration costs within each year 

τ_nb = copy(τ)

τ_nb[(N*C_τ*N+1):(2*N*C_τ*N),:] = τ[1:(N*C_τ*N),:]

ss_nb = steady_state(tol, maxit, N, C, R, ν[1, T], σ_0[1, T], σ_1[1, T], γ[1, T], B[:, T], s[:, T], τ_nb[:, T], r_bar[:, T], η[:, T], A[:, T], κ_0[:, T], κ_1[:, T], α[:, T], L_in[:, K], λ)

# the last period is the steady state
L_in_nb = copy(L_in)
L_in_nb[:, T] = ss_nb[:, 1]

V_in_nb = copy(V_in)
V_in_nb[:, T] = ss_nb[:, 2]

V_in_nb_2 = copy(V_in_2)
V_in_nb_2[:, T] = ss_nb[:, 2]

path_nb = transition_path(R, C, N, T, λ_2, L_in_nb, V_in_nb, V_in_nb_2, s, ν, τ_nb, M, σ_0, σ_1, κ_0, κ_1, A, r_bar, η, γ, B, tol, maxit)

V_nb = path_nb[:, 3:2+T]
L_nb = path_nb[:, (3+T):(2+2*T)]
real_wage_nb = path_nb[1:N*R*(C-1), (3+2*T):(2+3*T-1)]
nom_wage_nb =  path_nb[1:N*R*(C-1), (2+3*T):(1+4*T-1)]

gdp_nb = GDP(path_nb)

gdp_nb ./ gdp1

# common fertility across locations within races and ages 

# fertility per cohort-race-place 
α_data_c = CSV.read("../fertility/output/csv/alpha_common_matrix.csv")
α_mat_1_c = convert(Matrix, α_data_c)

# suppose that from period 9, age 20s (30) have the babies so that populations will be constant at the steady state 
# number of babies for nonblacks 
α_nb = 1/(s_mat2[1, 100] * s_mat2[2, 100])
α_b = 1/(s_mat2[8, 100] * s_mat2[9, 100])

# columns are periods.
# in rows:
# (i) races are in the upper tier,
# (ii) locations are in the middle tier,
# (iii) ages in the lower tier.

α_period_3 = vcat(repeat([0, 0, α_nb, 0, 0, 0, 0, 0], N), repeat([0, 0, α_b, 0, 0, 0, 0, 0], N))

α_mat_3 = repeat(α_period_3, 1, 100)

# 10 years of transition follow year 2010
# fertility in the transition period is between the fertility of 2010 and one of the steady state

α_mat_2_c = zeros(R*N*C, 10)

for i in 1:10
        α_mat_2_c[:, i] = (i/10) * α_period_3 + ((10-i)/10) * α_mat_1_c[:,K]
end
α_mat_2_c


α_c = hcat(α_mat_1_c, α_mat_2_c, α_mat_3)

ss_c = steady_state(tol, maxit, N, C, R, ν[1, T], σ_0[1, T], σ_1[1, T], γ[1, T], B[:, T], s[:, T], τ[:, T], r_bar[:, T], η[:, T], A[:, T], κ_0[:, T], κ_1[:, T], α_c[:, T], L_in[:, K], λ)

V_in = fill(1.0, (R*N*C), T)
V_in_2 = fill(2.0, (R*N*C), T)

# the last period is the steady state
L_in = copy(L_in)
L_in[:, T] = ss_c[:, 1]

V_in = copy(V_in)
V_in[:, T] = ss_c[:, 2]

V_in_2 = copy(V_in_2)
V_in_2[:, T] = ss_c[:, 2]

path_c = transition_path(R, C, N, T, λ_2, L_in, V_in, V_in_2, s, ν, τ, M, σ_0, σ_1, κ_0, κ_1, A, r_bar, η, γ, B, tol, maxit)

V_nb = path_nb[:, 3:2+T]
L_nb = path_nb[:, (3+T):(2+2*T)]
real_wage_nb = path_nb[1:N*R*(C-1), (3+2*T):(2+3*T-1)]
nom_wage_nb =  path_nb[1:N*R*(C-1), (2+3*T):(1+4*T-1)]

gdp_nb = GDP(path_nb)

gdp_nb ./ gdp1



############################################################################################
# heterogeneous rent elasticity

# place-specific shifter of rent 
r_bar_data2 = CSV.read("../rent_elas/output/csv/r_bar_2_hetero.csv")

r_bar_period2 = r_bar_data2.r_bar_2

r_bar2 = repeat(r_bar_period2, 1, T)

# place-specific elasticity of rent 
η_period2 = r_bar_data2.new_ivaw

η2 = repeat(η_period2, 1, T)

# steady state 

ss2 = steady_state(tol, maxit, N, C, R, ν[1, T], σ_0[1, T], σ_1[1, T], γ[1, T], B[:, T], s[:, T], τ[:, T], r_bar2[:, T], η2[:, T], A[:, T], κ_0[:, T], κ_1[:, T], α[:, T], L_in[:, K], λ)

V_mat_100_2 = reshape(ss2[:, 2], C, R*N)'
L_mat_100_2 = reshape(ss2[:, 1], C, R*N)'

real_wage_mat_100_2 = reshape(ss2[1:N*R*(C-1), 3], C-1, R*N)'

V_in = fill(1.0, (R*N*C), T)
V_in_2 = fill(2.0, (R*N*C), T)
L_in_2 = copy(L_in)

# the last period is the steady state
L_in_2[:, T] = ss2[:, 1]
V_in[:, T] = ss2[:, 2]
V_in_2[:, T] = ss2[:, 2]

path_h = transition_path(R, C, N, T, λ_2, L_in_2, V_in, V_in_2, s, ν, τ, M, σ_0, σ_1, κ_0, κ_1, A, r_bar2, η2, γ, B, tol, maxit)
Vh = path_h[:, 3:2+T]
Lh = path_h[:, (3+T):(2+2*T)]
real_wage_h = path_h[1:N*R*(C-1), (3+2*T):(2+3*T-1)]
nom_wage_h =  path_h[1:N*R*(C-1), (2+3*T):(1+4*T-1)]

gdph = GDP(path_h)

gdph ./ gdp1

# African Americans cannot move across the North and the South with hetero rent elasticity 
τ_nonorthsouth4060_b

ss_h_nonorthsouth4060_b = steady_state(tol, maxit, N, C, R, ν[1, T], σ_0[1, T], σ_1[1, T], γ[1, T], B[:, T], s[:, T], τ_nonorthsouth4060_b[:, T], r_bar2[:, T], η2[:, T], A[:, T], κ_0[:, T], κ_1[:, T], α[:, T], L_in[:, K], λ)

V_in = fill(1.0, (R*N*C), T)
V_in_2 = fill(2.0, (R*N*C), T)
L_in_2 = copy(L_in)

# the last period is the steady state
L_in_2[:, T] = ss_h_nonorthsouth4060_b[:, 1]
V_in[:, T] = ss_h_nonorthsouth4060_b[:, 2]
V_in_2[:, T] = ss_h_nonorthsouth4060_b[:, 2]

path_h_nonorthsouth4060_b = transition_path(R, C, N, T, λ_2, L_in_2, V_in, V_in_2, s, ν, τ_nonorthsouth4060_b, M, σ_0, σ_1, κ_0, κ_1, A, r_bar2, η2, γ, B, tol, maxit)
V_h_nonorthsouth4060_b = path_h_nonorthsouth4060_b[:, 3:2+T]
L_h_nonorthsouth4060_b = path_h_nonorthsouth4060_b[:, (3+T):(2+2*T)]
real_wage_h_nonorthsouth4060_b = path_h_nonorthsouth4060_b[1:N*R*(C-1), (3+2*T):(2+3*T-1)]
nom_wage_h_nonorthsouth4060_b =  path_h_nonorthsouth4060_b[1:N*R*(C-1), (2+3*T):(1+4*T-1)]

gdp_h_nonorthsouth4060_b = GDP(path_h_nonorthsouth4060_b)

gdp_h_nonorthsouth4060_b ./ gdph



# compute US income in data 
wage_data = CSV.read("../wage/state_race_cohort/output/csv/wage_q.csv")
w_mat_dt_1 = convert(Matrix, wage_data)
w_mat_dt_2 = repeat(w_mat_dt_1[:, K], 1, 100)
w_mat_dt = hcat(w_mat_dt_1, w_mat_dt_2)

(nom_wage1_big[:, 1:K] - w_mat_dt[:, 1:K]) ./ w_mat_dt[:, 1:K]
nom_wage1_big[:, 1:K]
w_mat_dt[:, 1:K]

income_dt = L_in .* w_mat_dt

US_income_dt = sum(income_dt, dims = 1)

w1_Ldt = L_in[:, 1:(T-1)] .* nom_wage1_big
US_w1_L_dt = sum(w1_Ldt, dims = 1)

wdt_L1 = w_mat_dt .* L1
US_wdt_L1 = sum(wdt_L1, dims = 1)



# double the migration costs from 1940 to 2010

τ_path3 = copy(τ)

τ_path3[:,1:8] = 2*τ[:,1:8]

τ_path3 

ss_path3 = steady_state(tol, maxit, N, C, R, ν[1, T], σ_0[1, T], σ_1[1, T], γ[1, T], B[:, T], s[:, T], τ_path3[:, T], r_bar[:, T], η[:, T], A[:, T], κ_0[:, T], κ_1[:, T], α[:, T], L_in[:, T], λ)

# the last period is the steady state
L_in_path3 = copy(L_in)
L_in_path3[:, T] = ss_path3[:, 1]

V_in_path3 = copy(V_in)
V_in_path3[:, T] = ss_path3[:, 2]

V_in_path3_2 = copy(V_in_2)
V_in_path3_2[:, T] = ss_path3[:, 2]

path3 = transition_path(R, C, N, T, λ_2, L_in_path3, V_in_path3, V_in_path3_2, s, ν, τ_path3, M, σ_0, σ_1, κ_0, κ_1, A, r_bar, η, γ, B, tol, maxit)

V3 = path3[:, 3:2+T]
L3 = path3[:, (3+T):(2+2*T)]
real_wage3 = path3[1:N*R*(C-1), (3+2*T):(2+3*T-1)]
nom_wage3 =  path3[1:N*R*(C-1), (2+3*T):(1+4*T-1)]

nom_wage3_big = zeros(R*N*C , (T-1))

for i in 1:(T-1)
        for j in 1:R 
                for n in 1:N 
                        nom_wage3_big[((j-1)*N*C+(n-1)*C+1):((j-1)*N*C+n*C) , i] = [0; nom_wage3[((j-1)*N*(C-1)+(n-1)*(C-1)+1):((j-1)*N*(C-1)+n*(C-1)), i]]
                end
        end
end

nom_wage3_big

income3 = L3[:, 1:(T-1)] .* nom_wage3_big

US_income3 = sum(income3, dims = 1)

years = 10 * (194:201)
plot(years, US_income1[1:K], legend=:topleft)
plot!(years, US_income_h[1:K], legend=:topleft)

plot!(years, US_w1_L_dt[1:K], legend=:topleft)
plot!(years, US_wdt_L1[1:K], legend=:topleft)
plot!(years, US_income2[1:K])
plot!(years, US_income3[1:K])
plot!(years, US_income_dt[1:K])

US_income2[1:K] ./ US_income1[1:K]

US_income3[1:K] ./ US_income1[1:K]

sum([1 3; 5 7], dims = 1)


nom_wage1
w_mat_dt

# compare populations from the model and those in data 
USpop1 = sum(L1, dims = 1)
USpop_dt = sum(L_in, dims = 1)

# compute the share of African Americans in the South 



L_in

south_pop_dt = zeros(2*R, T)
# columns are periods 
# rows are:
# (i) races are in the outer tier,
# (ii) north/south are in the inner tier

for t in 1:T 
        for i in 1:R 
                for n in 1:N 
                        northorsouth = north_or_south[n]
                        temp = L_in[((i-1)*N*C+(n-1)*C+1):((i-1)*N*C+n*C), t]
                        south_pop_dt[(i-1)*2+1+northorsouth , t] = south_pop_dt[(i-1)*2+1+northorsouth , t] + sum(temp)
                end
        end
end

south_pop_dt 

south_share_b_dt = south_pop_dt[4, :] ./ (south_pop_dt[3, :] + south_pop_dt[4, :])
south_share_nb_dt = south_pop_dt[2, :] ./ (south_pop_dt[1, :] + south_pop_dt[2, :])

years = 10 * (194:201)
plot(years, south_share_b[1:K], ylims = (0, 1))
plot!(years, south_share_b_h[1:K])
plot!(years, south_share_b_dt[1:K])


plot!(years, south_share_nb[1:K])
plot!(years, south_share_nb_h[1:K])
plot!(years, south_share_nb_dt[1:K])
