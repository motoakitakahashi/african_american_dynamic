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

# suppose the economy converges to a steady state in 100 periods from 2010
# from period 1 (1940) to period 8 (2010), I feed fundamentals backed out from data 
# from period 9, I fix fundamentals and the economy converges to the steady state 
# in period 108 (8+100), the economy reaches to the steady state

T = 108

# the number of places 
N = 38

# the number of cohorts in a given period
C = 8

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
amenity_mat2 = repeat(amenity_mat[:, 6], 1, 101) # suppose that amenities since 2010 are those in 2000
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

s_mat2 = repeat(s_mat[:, 8], 1, 100) # suppose that survival probabilities after 2010 are those in 2010

s = hcat(s_mat, s_mat2)
# columns are periods 
# each row is structured as:
# the survival probability for the youngest race 1, ..., the survival probability for the second oldest race 1, ...
# the survival probability for the youngest race 2, ..., the survival probability for the second oldest race 2, ...

# migration costs 
# note that the oldest can no longer migrate
τ_data = CSV.read("../mig_wage/output/csv/tau_matrix.csv")
# τ_data = CSV.read("../mig_wage/output/csv/tau_matrix_oridess.csv")
τ_mat = convert(Matrix, τ_data)

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

K_τ = 8
C_τ = 7

τ_1 = zeros(R*C_τ*N*N, K_τ)

for i in 1:K_τ
        temp = τ_mat[((i-1)R*N+1):(i*R*N) , :]'
        temp = reshape(temp, R*C_τ*N*N)
        τ_1[:, i] = temp

end

τ_2 = repeat(τ_1[:, K_τ], 1, 100)

τ = hcat(τ_1, τ_2)


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

A_mat_2 = repeat(A_mat_1[:, K], 1, 100)


A = hcat(A_mat_1, A_mat_2)

A[:, 2] - A[:, 1]

# cohort-specific productivity 
κ_0_data = CSV.read("../elas_sub/output/csv/age_prod_matrix.csv")
κ_0_mat_1 = convert(Matrix, κ_0_data)

κ_0_mat_2 = repeat(κ_0_mat_1[:, K], 1, 100)


κ_0 = hcat(κ_0_mat_1, κ_0_mat_2)

# κ_0_period = [1/7, 1/7, 1/7, 1/7, 1/7, 1/7, 1/7,
#               1/7, 1/7, 1/7, 1/7, 1/7, 1/7, 1/7]
# the relative productivity of the second youngest in place 1, ..., the relative productivity of the oldest in place 1,
# the relative productivity of the second youngest in place 2, ..., the relative productivity of the oldest in place 2



# race-specific productivity (within cohorts)
κ_1_data = CSV.read("../elas_sub/output/csv/race_prod_matrix.csv")
κ_1_mat_1 = convert(Matrix, κ_1_data)

κ_1_mat_2 = repeat(κ_1_mat_1[:, K], 1, 100)

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

α_period_2 = vcat(repeat([0, 0, α_nb, 0, 0, 0, 0, 0], N), repeat([0, 0, α_b, 0, 0, 0, 0, 0], N))

α_mat_2 = repeat(α_period_2, 1, 100)

α = hcat(α_mat_1, α_mat_2)

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

M = repeat(zeros(R*C*N) , 1, 108)

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
L_mat_2 = repeat(L_mat_1[:, K], 1, 100)

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

ss = steady_state(tol, maxit, N, C, R, ν[1, T], σ_0[1, T], σ_1[1, T], γ[1, T], B[:, T], s[:, T], τ[:, T], r_bar[:, T], η[:, T], A[:, T], κ_0[:, T], κ_1[:, T], α[:, T], L_in[:, T], λ)

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


function transition_path(R, C, N, T, λ_2, L_in, V_in, V_in_2, s, ν, τ, M, σ_0, σ_1, κ_0, κ_1, A, r_bar, η, γ, B, tol, maxit)
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

# compute GDPs
nom_wage1_big = zeros(R*N*C , (T-1))

for i in 1:(T-1)
        for j in 1:R 
                for n in 1:N 
                        nom_wage1_big[((j-1)*N*C+(n-1)*C+1):((j-1)*N*C+n*C) , i] = [0; nom_wage1[((j-1)*N*(C-1)+(n-1)*(C-1)+1):((j-1)*N*(C-1)+n*(C-1)), i]]
                end
        end
end

nom_wage1_big

income1 = L1[:, 1:(T-1)] .* nom_wage1_big

US_income1 = sum(income1, dims = 1)

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

ss2 = steady_state(tol, maxit, N, C, R, ν[1, T], σ_0[1, T], σ_1[1, T], γ[1, T], B[:, T], s[:, T], τ[:, T], r_bar2[:, T], η2[:, T], A[:, T], κ_0[:, T], κ_1[:, T], α[:, T], L_in[:, T], λ)

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

# compute GDPs
nom_wage_h_big = zeros(R*N*C , (T-1))

for i in 1:(T-1)
        for j in 1:R 
                for n in 1:N 
                        nom_wage_h_big[((j-1)*N*C+(n-1)*C+1):((j-1)*N*C+n*C) , i] = [0; nom_wage_h[((j-1)*N*(C-1)+(n-1)*(C-1)+1):((j-1)*N*(C-1)+n*(C-1)), i]]
                end
        end
end

nom_wage_h_big

income_h = Lh[:, 1:(T-1)] .* nom_wage_h_big

US_income_h = sum(income_h, dims = 1)



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

# counterfactual economy in which migration costs are zero 
τ_path2 = zeros(R*C_τ*N*N, 108)

ss_path2 = steady_state(tol, maxit, N, C, R, ν[1, T], σ_0[1, T], σ_1[1, T], γ[1, T], B[:, T], s[:, T], τ_path2[:, T], r_bar[:, T], η[:, T], A[:, T], κ_0[:, T], κ_1[:, T], α[:, T], L_in[:, T], λ)

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

nom_wage2_big = zeros(R*N*C , (T-1))

for i in 1:(T-1)
        for j in 1:R 
                for n in 1:N 
                        nom_wage2_big[((j-1)*N*C+(n-1)*C+1):((j-1)*N*C+n*C) , i] = [0; nom_wage2[((j-1)*N*(C-1)+(n-1)*(C-1)+1):((j-1)*N*(C-1)+n*(C-1)), i]]
                end
        end
end

nom_wage2_big

income2 = L2[:, 1:(T-1)] .* nom_wage2_big

US_income2 = sum(income2, dims = 1)

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

# read the concordance table mapping geographic units to the South 

geo_south = CSV.read("../geo_unit/output/my_geo_south.csv")

south_ind = geo_south.south
# takes 0 if the location is in the North,
# takes 1 if the location is in the South.

south_pop = zeros(2*R, T)
south_pop_h = zeros(2*R, T)
# columns are periods 
# rows are:
# (i) races are in the outer tier,
# (ii) north/south are in the inner tier

for t in 1:T 
        for i in 1:R 
                for n in 1:N 
                        northorsouth = south_ind[n]
                        temp = L1[((i-1)*N*C+(n-1)*C+1):((i-1)*N*C+n*C), t]
                        temp2 = Lh[((i-1)*N*C+(n-1)*C+1):((i-1)*N*C+n*C), t]
                        south_pop[(i-1)*2+1+northorsouth , t] = south_pop[(i-1)*2+1+northorsouth , t] + sum(temp)
                        south_pop_h[(i-1)*2+1+northorsouth , t] = south_pop_h[(i-1)*2+1+northorsouth , t] + sum(temp2)
                end
        end
end

south_share_b = south_pop[4, :] ./ (south_pop[3, :] + south_pop[4, :])
south_share_nb = south_pop[2, :] ./ (south_pop[1, :] + south_pop[2, :])

south_share_b_h = south_pop_h[4, :] ./ (south_pop_h[3, :] + south_pop_h[4, :])
south_share_nb_h = south_pop_h[2, :] ./ (south_pop_h[1, :] + south_pop_h[2, :])

L_in

south_pop_dt = zeros(2*R, T)
# columns are periods 
# rows are:
# (i) races are in the outer tier,
# (ii) north/south are in the inner tier

for t in 1:T 
        for i in 1:R 
                for n in 1:N 
                        northorsouth = south_ind[n]
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
