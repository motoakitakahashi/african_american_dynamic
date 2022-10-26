# Motoaki Takahashi 

# Julia 1.7.2

# March, April, August, September, October 2022

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
# cd("D:/onedrive/OneDrive - The Pennsylvania State University/dynamic/code")

# for the laptop
cd("C:/Users/takah/OneDrive - The Pennsylvania State University/dynamic/code")

# suppose the economy converges to a steady state in 100 periods from 2010
# from period 1 (1940) to period 8 (2010), I feed fundamentals backed out from data 
# from period 9, I fix fundamentals and the economy converges to the steady state 
# in period 108 (8+100), the economy reaches the steady state

T = 108

# the number of places 
N = 38

# years in the sample 
years = 1940:10:2010

# the number of cohorts in a given period
C = 8

# the number of races 
R = 2

# the number of years in data (1940, 1950, ..., 2010)
K = 8

# for parameters, columns generically represent periods,
# rows are races/places 

# migration elesticity
# ν = fill(1/1.499, 1, T)
ν = fill(1/1.38428  , 1, T)

# Cobb-Douglas share of housing
γ = fill(0.25, 1, T)

# amenities 
amenity_data = CSV.read("../mig_wage/output/csv/amenity_matrix_stay.csv")
amenity_data_d = CSV.read("../mig_wage/output/csv/amenity_matrix_d_stay.csv")


amenity_mat = convert(Matrix, amenity_data)
amenity_mat_d = convert(Matrix, amenity_data_d)

# amenity_1940 = amenity_mat[:, 1] # suppose that amenities in 1940 are the same as those in 1950
amenity_1940 = amenity_mat_d[:, 1]
amenity_mat2 = repeat(amenity_mat[:, 7], 1, (T-8)) # suppose that amenities since 2010 are those in 2000
B = hcat(amenity_1940, amenity_mat, amenity_mat2)

B_d = hcat(amenity_mat_d, repeat(amenity_mat_d[:, 8], 1, (T-8)))

# columns are periods
# each row is structured as follows:
# the amenity for the youngest race 1 in place 1, ..., the amenity for the oldest race 1 in place 1, ...
# the amenity for the youngest race 1 in place 2, ..., the amenity for the oldest race 1 in place 2, ...
# the amenity for the youngest race 2 in place 1, ..., the amenity for the oldest race 2 in place 1, ...
# the amenity for the youngest race 2 in place 2, ..., the amenity for the oldest race 2 in place 2, ...



# survival probabilities
s_data = CSV.read("../life_table/output/s_matrix_q.csv")
s_mat = convert(Matrix, s_data)

# remove survival probabilities for the ages 70s 
# s_mat = s_mat[[1:(C-1);(C+1):(2*C-1)], :]


s_mat2 = repeat(s_mat[:, 8], 1, (T-K)) # suppose that survival probabilities after 2010 are those in 2010

s = hcat(s_mat, s_mat2)
# columns are periods 
# each row is structured as:
# the survival probability for the youngest race 1, ..., the survival probability for the second oldest race 1, ...
# the survival probability for the youngest race 2, ..., the survival probability for the second oldest race 2, ...

# migration costs 
# note that the oldest can no longer migrate
# τ_data = CSV.read("../mig_wage/output/csv/tau_matrix_sf2.csv")
# τ_data = CSV.read("../mig_wage/output/csv/tau_matrix.csv")
# τ_data = CSV.read("../mig_wage/output/csv/tau_matrix_oridess.csv")
τ_data = CSV.read("../mig_wage/output/csv/tau_matrix_stay.csv")

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

# fertility per cohort-race-place 
α_data = CSV.read("../fertility/output/csv/alpha_common_matrix.csv")

# α_data = CSV.read("../fertility/output/csv/alpha_matrix.csv")

α_mat_1 = convert(Matrix, α_data)

# suppose that from period 9, age 20s (30) have the babies so that populations will be constant at the steady state 
# number of babies for nonblacks 
α_nb = 1/(s_mat2[1, 100] * s_mat2[2, 100])
α_b = 1/(s_mat2[8, 100] * s_mat2[9, 100])

# α_nb = 1/s_mat2[1, 100]
# α_b = 1/s_mat2[7, 100]

# columns are periods.
# in rows:
# (i) races are in the upper tier,
# (ii) locations are in the middle tier,
# (iii) ages in the lower tier.

α_period_2 = vcat(repeat([0, 0, α_nb, 0, 0, 0, 0, 0], N), repeat([0, 0, α_b, 0, 0, 0, 0, 0], N))

α_mat_2 = repeat(α_period_2, 1, (T-K+1))

# from 2010, fertility rate is such that populations are constant

α = hcat(α_mat_1[:, 1:(K-1)], α_mat_2)

# I guess I need to restrict the value of α to get a steady state, 
# α_period = [    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
#                 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
#                 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
#                 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,]
# the fertility of the second youngest race 1 in place 1, ..., the fertility of the oldest race 1 in place 1,
# the fertility of the second youngest race 1 in place 2, ..., the fertility of the oldest race 1 in place 2, ???
# the fertility of the second youngest race 2 in place 1, ..., the fertility of the oldest race 2 in place 1,
# the fertility of the second youngest race 2 in place 2, ..., the fertility of the oldest race 2 in place 2



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
            .* repeat(L_cohort_mat, R, 1) .^ (-1/σ_0 + 1/σ_1) .* κ_1_mat .^ (1/σ_1) .* L_mat[:, 2:C] .^ (-1/σ_1))


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

L_in_mat = reshape(L_in[:, K], C, (R*N))'
s[:, K]
sum(L_in_mat[1:N, 2])
sum(L_in_mat[1:N, 3])

sum(L_in_mat[(N+1):R*N, 2])
sum(L_in_mat[(N+1):R*N, 3])

sum(L_in_mat[1:N, 3]) / sum(L_in_mat[1:N, 2])

sum(L_in_mat[(N+1):R*N, 3]) / sum(L_in_mat[(N+1):R*N, 2])

# the initial period (0) is year 1940 from data

# rows are race-place-age-in-the-current-periods, columns are periods
# in each column (period), 
# the 1st row is the population of the youngest race 1 in place 1,
# the 2nd row is the population of the second youngest race 1 in place 1,
# ...,
# the 32nd row is the population of the oldest race 2 in place 2

tol = 10 ^ (-7)
maxit = 5000
tol_2 = 10 ^ (-6)

# the dumpening parameter for iteration
λ = 0.5

function steady_state(tol, maxit, N, C, R, ν, σ_0, σ_1, γ, B, s, τ, r_bar, η, A, κ_0, κ_1, α, L_in, λ)
        # L_in_mat = reshape(L_in, C, (R*N))'
        # sum_o_30 = sum(L_in_mat[1:N, 3])
        # sum_b_30 = sum(L_in_mat[(N+1):R*N, 3])


        L = copy(L_in)
        V = []
        w = []
        r = []
        real_wage = []

        dif = 1.0
        count = 0
        # V_mat = zeros((R*N), C)
        # eventually output = reshape(V_mat', (R*N*C) )
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
                println("dif: ", dif)

                L = (1-λ) * L + λ * L_new

                # L_mat = reshape(L, C, (R*N))'
        
                # multiplier_o = sum_o_30 / sum(L_mat[1:N, 3])
                # multiplier_b = sum_b_30 / sum(L_mat[(N+1):R*N, 3])

                # L_mat[1:N, 3] = multiplier_o * L_mat[1:N, 3]
                # L_mat[(N+1):R*N, 3] = multiplier_b * L_mat[(N+1):R*N, 3]

                # L = reshape(L_mat', (R*N*C))


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

# I compute populations forward, given expected values.
# Then given populations, I update expected values.
λ_2 = 0.5

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
                println("dif: ", dif)
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


function ss_tp(R, C, N, T, λ, λ_2, L_in, s, ν, τ, M, σ_0, σ_1, κ_0, κ_1, A, r_bar, η, γ, B, tol, tol_2, maxit, α)
        ss = steady_state(tol, maxit, N, C, R, ν[1, T], σ_0[1, T], σ_1[1, T], γ[1, T], B[:, T], s[:, T], τ[:, T], r_bar[:, T], η[:, T], A[:, T], κ_0[:, T], κ_1[:, T], α[:, T], L_in[:, K], λ)

        # arbitrary initial guesses for expected values
        V_in = fill(1.0, (R*N*C), T)
        V_in_2 = fill(2.0, (R*N*C), T)

        # the last period is the steady state
        L_in[:, T] = ss[:, 1]
        V_in[:, T] = ss[:, 2]
        V_in_2[:, T] = ss[:, 2]

        path = transition_path(R, C, N, T, λ_2, L_in, V_in, V_in_2, s, ν, τ, M, σ_0, σ_1, κ_0, κ_1, A, r_bar, η, γ, B, tol_2, maxit, α)

        return path
end

# the following function maps a transition path of nominal wages and populations to 
# nation-wide GDPs and per capita GDPs
function GDP(path_res, T)
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

function local_pop(L, T, R, N, C)
        output = zeros(R*N, T)

        for i in 1:R 
                for n in 1:N 
                        output[(i-1)*N+n, :] = sum(L[((i-1)*N*C+(n-1)*C+1):((i-1)*N*C+n*C), :], dims = 1)
                end
        end
        return output
end

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

union_or_conf = south_dt.conf

conf_ind_df = filter(row -> row.conf == 1, south_dt)
conf_ind = conf_ind_df.x1

union_ind_df = filter(row -> row.conf == 0, south_dt)
union_ind = union_ind_df.x1

N_C = length(conf_ind)
N_U = length(union_ind)

function share_south(L, north_or_south, T)

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

# wage matrix including the youngest (who are not working)

function wage_big(wage_small)
        output = zeros(R*N*C , (T-1))

        for i in 1:(T-1)
                for j in 1:R 
                        for n in 1:N 
                                output[((j-1)*N*C+(n-1)*C+1):((j-1)*N*C+n*C) , i] = [0; wage_small[((j-1)*N*(C-1)+(n-1)*(C-1)+1):((j-1)*N*(C-1)+n*(C-1)), i]]
                        end
                end
        end
        return output
end

# compute US GDPs in data
wage_data = CSV.read("../wage/state_race_cohort/output/csv/wage_q.csv")
w_mat_dt_1 = convert(Matrix, wage_data)
w_mat_dt_2 = repeat(w_mat_dt_1[:, K], 1, 100)
w_mat_dt = hcat(w_mat_dt_1[:, 1:K], w_mat_dt_2)
income_dt = L_in .* w_mat_dt

US_income_dt = sum(income_dt, dims = 1)

per_capita_US_income_dt = US_income_dt ./ sum(L_in, dims = 1)

local_pop_dt = local_pop(L_in, T, R, N, C)

share_south_dt = share_south(L_in, north_or_south, T)

share_conf_dt = share_south(L_in, union_or_conf, T)


# LEVEL ESTIMATION

# IV 1

# elasticity of substitution between ages 
σ_0_level = fill(1/0.4322, 1, T) 
# 0.3401 without 80
# 0.4322 with 80


# elasticity of substitution between races (within ages)
σ_1_level = fill(1/0.1184, 1, T) 
# 0.1108 without 80
# 0.1184 with 80
# productivity
A_data_level = CSV.read("../elas_sub/output/csv/loc_prod_matrix_level.csv")

A_mat_1_level =  convert(Matrix, A_data_level)

A_mat_2_level = repeat(A_mat_1_level[:, K], 1, (T-K))

A_level = hcat(A_mat_1_level, A_mat_2_level)

# cohort-specific productivity 
κ_0_data_level = CSV.read("../elas_sub/output/csv/age_prod_matrix_level.csv")
κ_0_mat_1_level = convert(Matrix, κ_0_data_level)

κ_0_mat_2_level = repeat(κ_0_mat_1_level[:, K], 1, (T-K))

κ_0_level = hcat(κ_0_mat_1_level, κ_0_mat_2_level)

# race-specific productivity (within cohorts)
κ_1_data_level = CSV.read("../elas_sub/output/csv/race_prod_matrix_level.csv")
κ_1_mat_1_level = convert(Matrix, κ_1_data_level)

κ_1_mat_2_level = repeat(κ_1_mat_1_level[:, K], 1, (T-K))

κ_1_level = hcat(κ_1_mat_1_level, κ_1_mat_2_level)

# immgrants from abroad 
M_data = CSV.read("../migration_flow/state/output/csv/immigrant.csv") # this data contains the numbers of immigrants from 1950 to 2010
M_mat_1 = convert(Matrix, M_data)

# I assume that since period 9 (2020), noone comes to the US from abroad, and noone leaves the US to abroad

M_mat_2 = repeat(zeros(R*C*N) , 1, 100)

M1 = hcat(zeros(R*C*N), M_mat_1, M_mat_2)

# alternatively, I can assume zero immigrants from abroad
M0 = repeat(zeros(R*C*N) , 1, T)

# baseline equilibrium with immigrants 
path1_M1_level = ss_tp(R, C, N, T, λ, λ_2, L_in, s, ν, τ, M1, σ_0_level, σ_1_level, κ_0_level, κ_1_level, A_level, r_bar, η, γ, B, tol, tol_2, maxit, α)


V1_M1_level = path1_M1_level[:, 3:2+T]


L1_M1_level = path1_M1_level[:, (3+T):(2+2*T)]
nom_wage1_M1_level =  path1_M1_level[1:N*R*(C-1), (2+3*T):(1+4*T-1)]

gdp1_M1_level = GDP(path1_M1_level, T)
local_gdp1_M1_level = local_GDP(path1_M1_level, N, T)

share_south1_M1_level = share_south(L1_M1_level, north_or_south, T)
share_conf_M1_level = share_south(L1_M1_level, union_or_conf, T)

local_pop_M1_level = local_pop(L1_M1_level, T, R, N, C)

V_youngest1_M1_level = get_V_youngest(V1_M1_level)

# from 1940 to 1960, only African Americans have choke migration costs between the Union and the confederacy
τ_nounionconf4060_b_data = CSV.read("../mig_wage/output/csv/tau_nounionconf4060_b_matrix_stay.csv")

τ_nounionconf4060_b = make_τ_mat(τ_nounionconf4060_b_data, K_τ, C_τ)

path_nounionconf4060_b_M1_level = ss_tp(R, C, N, T, λ, λ_2, L_in, s, ν, τ_nounionconf4060_b, M1, σ_0_level, σ_1_level, κ_0_level, κ_1_level, A_level, r_bar, η, γ, B, tol, tol_2, maxit, α)

V_nounionconf4060_b_M1_level = path_nounionconf4060_b_M1_level[:, 3:2+T]

gdp_nounionconf4060_b_M1_level = GDP(path_nounionconf4060_b_M1_level, T)

gdp_ratio_nounionconf4060_b_M1_level = gdp_nounionconf4060_b_M1_level ./ gdp1_M1_level

gdp1_M1_level ./ gdp_nounionconf4060_b_M1_level

V_nounionconf4060_b_M1_level_youngest = get_V_youngest(V_nounionconf4060_b_M1_level)
V_nounionconf4060_b_M1_level_youngest_ratio = V_nounionconf4060_b_M1_level_youngest ./ V_youngest1_M1_level

save_V_nounionconf4060_b_M1_level_youngest_ratio = hcat(V_nounionconf4060_b_M1_level_youngest_ratio[1:N, 1],
V_nounionconf4060_b_M1_level_youngest_ratio[(N+1):R*N, 1])

save_V_nounionconf4060_b_M1_level_youngest_ratio = convert(DataFrame, save_V_nounionconf4060_b_M1_level_youngest_ratio)
CSV.write("output/csv/V_nounionconf4060_b_M1_level_youngest_ratio_stay.csv", save_V_nounionconf4060_b_M1_level_youngest_ratio)

save_gdp_ratio_nounionconf4060_b_M1_level = convert(DataFrame, gdp_ratio_nounionconf4060_b_M1_level)
CSV.write("output/csv/gdp_ratio_nounionconf4060_b_M1_level_stay.csv", save_gdp_ratio_nounionconf4060_b_M1_level)


# from 1940 to 1960, only nonblacks have choke migration costs between the union and the Confederacy
τ_nounionconf4060_o_data = CSV.read("../mig_wage/output/csv/tau_nounionconf4060_o_matrix_stay.csv")

τ_nounionconf4060_o = make_τ_mat(τ_nounionconf4060_o_data, K_τ, C_τ)

path_nounionconf4060_o_M1_level = ss_tp(R, C, N, T, λ, λ_2, L_in, s, ν, τ_nounionconf4060_o, M1, σ_0_level, σ_1_level, κ_0_level, κ_1_level, A_level, r_bar, η, γ, B, tol, tol_2, maxit, α)

V_nounionconf4060_o_M1_level = path_nounionconf4060_o_M1_level[:, 3:2+T]

gdp_nounionconf4060_o_M1_level = GDP(path_nounionconf4060_o_M1_level, T)

gdp_ratio_nounionconf4060_o_M1_level = gdp_nounionconf4060_o_M1_level ./ gdp1_M1_level

V_nounionconf4060_o_M1_level_youngest = get_V_youngest(V_nounionconf4060_o_M1_level)
V_nounionconf4060_o_M1_level_youngest_ratio = V_nounionconf4060_o_M1_level_youngest ./ V_youngest1_M1_level

save_V_nounionconf4060_o_M1_level_youngest_ratio = hcat(V_nounionconf4060_o_M1_level_youngest_ratio[1:N, 1],
V_nounionconf4060_o_M1_level_youngest_ratio[(N+1):R*N, 1])

save_V_nounionconf4060_o_M1_level_youngest_ratio = convert(DataFrame, save_V_nounionconf4060_o_M1_level_youngest_ratio)
CSV.write("output/csv/V_nounionconf4060_o_M1_level_youngest_ratio_stay.csv", save_V_nounionconf4060_o_M1_level_youngest_ratio)


# No migration from the confederacy to the union for African Americans (suggested by Jonathan)
τ_no_conf_to_union_4060_b_data = CSV.read("../mig_wage/output/csv/tau_no_conf_to_union_4060_b_matrix_stay.csv")

τ_no_conf_to_union_4060_b = make_τ_mat(τ_no_conf_to_union_4060_b_data, K_τ, C_τ)

path_no_conf_to_union_4060_b_M1_level = ss_tp(R, C, N, T, λ, λ_2, L_in, s, ν, τ_no_conf_to_union_4060_b, M1, σ_0_level, σ_1_level, κ_0_level, κ_1_level, A_level, r_bar, η, γ, B, tol, tol_2, maxit, α)

V_no_conf_to_union_4060_b_M1_level = path_no_conf_to_union_4060_b_M1_level[:, 3:2+T]

gdp_no_conf_to_union_4060_b_M1_level = GDP(path_no_conf_to_union_4060_b_M1_level, T)

gdp_ratio_no_conf_to_union_4060_b_M1_level = gdp_no_conf_to_union_4060_b_M1_level ./ gdp1_M1_level

V_no_conf_to_union_4060_b_M1_level_youngest = get_V_youngest(V_no_conf_to_union_4060_b_M1_level)
V_no_conf_to_union_4060_b_M1_level_youngest_ratio = V_no_conf_to_union_4060_b_M1_level_youngest ./ V_youngest1_M1_level

save_V_no_conf_to_union_4060_b_M1_level_youngest_ratio = hcat(V_no_conf_to_union_4060_b_M1_level_youngest_ratio[1:N, 1],
V_no_conf_to_union_4060_b_M1_level_youngest_ratio[(N+1):R*N, 1])

save_V_no_conf_to_union_4060_b_M1_level_youngest_ratio = convert(DataFrame, save_V_no_conf_to_union_4060_b_M1_level_youngest_ratio)
CSV.write("output/csv/V_no_conf_to_union_4060_b_M1_level_youngest_ratio_stay.csv", save_V_no_conf_to_union_4060_b_M1_level_youngest_ratio)


# No migration from the confederacy to the union for the others (suggested by Jonathan)
τ_no_conf_to_union_4060_o_data = CSV.read("../mig_wage/output/csv/tau_no_conf_to_union_4060_o_matrix_stay.csv")

τ_no_conf_to_union_4060_o = make_τ_mat(τ_no_conf_to_union_4060_o_data, K_τ, C_τ)

path_no_conf_to_union_4060_o_M1_level = ss_tp(R, C, N, T, λ, λ_2, L_in, s, ν, τ_no_conf_to_union_4060_o, M1, σ_0_level, σ_1_level, κ_0_level, κ_1_level, A_level, r_bar, η, γ, B, tol, tol_2, maxit, α)

V_no_conf_to_union_4060_o_M1_level = path_no_conf_to_union_4060_o_M1_level[:, 3:2+T]

gdp_no_conf_to_union_4060_o_M1_level = GDP(path_no_conf_to_union_4060_o_M1_level, T)

gdp_ratio_no_conf_to_union_4060_o_M1_level = gdp_no_conf_to_union_4060_o_M1_level ./ gdp1_M1_level

V_no_conf_to_union_4060_o_M1_level_youngest = get_V_youngest(V_no_conf_to_union_4060_o_M1_level)
V_no_conf_to_union_4060_o_M1_level_youngest_ratio = V_no_conf_to_union_4060_o_M1_level_youngest ./ V_youngest1_M1_level

save_V_no_conf_to_union_4060_o_M1_level_youngest_ratio = hcat(V_no_conf_to_union_4060_o_M1_level_youngest_ratio[1:N, 1],
V_no_conf_to_union_4060_o_M1_level_youngest_ratio[(N+1):R*N, 1])

save_V_no_conf_to_union_4060_o_M1_level_youngest_ratio = convert(DataFrame, save_V_no_conf_to_union_4060_o_M1_level_youngest_ratio)
CSV.write("output/csv/V_no_conf_to_union_4060_o_M1_level_youngest_ratio_stay.csv", save_V_no_conf_to_union_4060_o_M1_level_youngest_ratio)


# No migration from the Union to the Confederacy for the others (suggested by Jonathan)

τ_no_union_to_conf_4060_o_data = CSV.read("../mig_wage/output/csv/tau_no_union_to_conf_4060_o_matrix_stay.csv")

τ_no_union_to_conf_4060_o = make_τ_mat(τ_no_union_to_conf_4060_o_data, K_τ, C_τ)

path_no_union_to_conf_4060_o_M1_level = ss_tp(R, C, N, T, λ, λ_2, L_in, s, ν, τ_no_union_to_conf_4060_o, M1, σ_0_level, σ_1_level, κ_0_level, κ_1_level, A_level, r_bar, η, γ, B, tol, tol_2, maxit, α)

V_no_union_to_conf_4060_o_M1_level = path_no_union_to_conf_4060_o_M1_level[:, 3:2+T]

gdp_no_union_to_conf_4060_o_M1_level = GDP(path_no_union_to_conf_4060_o_M1_level, T)

gdp_ratio_no_union_to_conf_4060_o_M1_level = gdp_no_union_to_conf_4060_o_M1_level ./ gdp1_M1_level

V_no_union_to_conf_4060_o_M1_level_youngest = get_V_youngest(V_no_union_to_conf_4060_o_M1_level)
V_no_union_to_conf_4060_o_M1_level_youngest_ratio = V_no_union_to_conf_4060_o_M1_level_youngest ./ V_youngest1_M1_level

save_V_no_union_to_conf_4060_o_M1_level_youngest_ratio = hcat(V_no_union_to_conf_4060_o_M1_level_youngest_ratio[1:N, 1],
V_no_union_to_conf_4060_o_M1_level_youngest_ratio[(N+1):R*N, 1])

save_V_no_union_to_conf_4060_o_M1_level_youngest_ratio = convert(DataFrame, save_V_no_union_to_conf_4060_o_M1_level_youngest_ratio)
CSV.write("output/csv/V_no_union_to_conf_4060_o_M1_level_youngest_ratio_stay.csv", save_V_no_union_to_conf_4060_o_M1_level_youngest_ratio)


# draw graphs 

# compare US labor incomes 
plot(years, US_income_dt[1, 1:K], legend = :topleft, label = "data", c = :black,
xlab = "year", ylab = "US labor income")
plot!(years, gdp1_M1_level[1, 1:K], label = "model", c=:black, line=:dash)
# savefig("output/figures/model_fit/income_model_data_stay.pdf")
gdp1_M1_level[1, 1:K] ./ US_income_dt[1, 1:K]

# compare population vectors at location-age-race levels
pop_cor = zeros(K)
for i in 1:K
        pop_cor[i] = cor(L_in[:, i], L1_M1_level[:, i])
end

plot(years, pop_cor, legend = false, xlab = "year", ylab = "correlation",
c = :black, ylims = [0.90, 1.02])
# savefig("output/figures/model_fit/correlation_pop_model_data_stay.pdf")

# compare populations vectors at location-race lavels 
pop_cor_loc = zeros(K)
for i in 1:K
        pop_cor_loc[i] = cor(local_pop_dt[:, i], local_pop_M1_level[:, i])
end

plot(years, pop_cor_loc, legend = false, xlab = "year", ylab = "correlation",
c = :black, ylims = [0.90, 1.02])
# savefig("output/figures/model_fit/correlation_pop_model_data_loc_stay.pdf")

plot(years, gdp_ratio_no_conf_to_union_4060_b_M1_level[1, 1:K], xlabel = "year",
        ylabel ="US income relative to the baseline", label = "No migraiton from the South to the North", c=:black,
        ylims = [0.989, 1.001],
        legend = :topright)
plot!(years, gdp_ratio_nounionconf4060_b_M1_level[1, 1:K],
 label = "No migraiton across the North and the South", c=:black,
        line = :dash)
# savefig("output/figures/counterfactual/income_no_mig_4060_b_union_conf_bilateral_vs_unilateral_stay.pdf")

plot(years - fill(10.0, K), V_nounionconf4060_b_M1_level_youngest_ratio[N+21, 1:K], xlabel = "cohort",
ylabel = "initial expected value", label = "Mississippi", legend = :bottomright, c=:black)
plot!(years - fill(10.0, K), V_nounionconf4060_b_M1_level_youngest_ratio[N+10, 1:K], label = "Georgia", c=:black, line = :dash)
plot!(years - fill(10.0, K), V_nounionconf4060_b_M1_level_youngest_ratio[N+11, 1:K], label = "Illinois", c=:blue, line = :solid)
plot!(years - fill(10.0, K), V_nounionconf4060_b_M1_level_youngest_ratio[N+25, 1:K], label = "New York", c=:blue, line = :dash)
# savefig("output/figures/counterfactual/expected_value_cohort_V_nounionconf4060_b_stay.pdf")

plot(years - fill(10.0, K), V_nounionconf4060_b_M1_level_youngest_ratio[N+21, 1:K], xlabel = "cohort",
ylabel = "initial expected value", label = "Mississippi: African Americans", legend = :bottomright, c=:black)
plot!(years - fill(10.0, K), V_nounionconf4060_b_M1_level_youngest_ratio[21, 1:K], label = "Mississippi: Others", c=:black, line = :dash)
plot!(years - fill(10.0, K), V_nounionconf4060_b_M1_level_youngest_ratio[N+11, 1:K], label = "Illinois: African Americans", c=:blue, line = :solid)
plot!(years - fill(10.0, K), V_nounionconf4060_b_M1_level_youngest_ratio[11, 1:K], label = "Illinois: Others", c=:blue, line = :dash)
# savefig("output/figures/counterfactual/expected_value_cohort_V_nounionconf4060_b_bo_stay.pdf")

plot(years - fill(10.0, K), V_no_conf_to_union_4060_b_M1_level_youngest_ratio[N+21, 1:K], xlabel = "cohort",
ylabel = "initial expected value", label = "Mississippi", legend = :bottomright, c=:black)
plot!(years - fill(10.0, K), V_no_conf_to_union_4060_b_M1_level_youngest_ratio[N+10, 1:K], label = "Georgia", c=:black, line = :dash)
plot!(years - fill(10.0, K), V_no_conf_to_union_4060_b_M1_level_youngest_ratio[N+11, 1:K], label = "Illinois", c=:blue, line = :solid)
plot!(years - fill(10.0, K), V_no_conf_to_union_4060_b_M1_level_youngest_ratio[N+25, 1:K], label = "New York", c=:blue, line = :dash)
savefig("output/figures/counterfactual/expected_value_cohort_V_no_conf_to_union_4060_b_stay.pdf")

plot(years - fill(10.0, K), V_no_conf_to_union_4060_b_M1_level_youngest_ratio[N+21, 1:K], xlabel = "cohort",
ylabel = "initial expected value", label = "Mississippi: African Americans", legend = :bottomright, c=:black)
plot!(years - fill(10.0, K), V_no_conf_to_union_4060_b_M1_level_youngest_ratio[21, 1:K],
 label = "Mississippi: Others", c=:black, line = :dash)
plot!(years - fill(10.0, K), V_no_conf_to_union_4060_b_M1_level_youngest_ratio[N+11, 1:K],
label = "Illinois: African Americans", c=:blue, line = :solid)
plot!(years - fill(10.0, K), V_no_conf_to_union_4060_b_M1_level_youngest_ratio[11, 1:K],
label = "Illinois: Others", c=:blue, line = :dash)
savefig("output/figures/counterfactual/expected_value_cohort_V_no_conf_to_union_4060_b_bo_stay.pdf")

plot(years, gdp_ratio_no_union_to_conf_4060_o_M1_level[1, 1:K], xlabel = "year",
        ylabel ="US income relative to the baseline", label = "No migraiton from the North to the South", c=:black,
        ylims = [0.99, 1.01],
        legend = :bottomright)
plot!(years, gdp_ratio_nounionconf4060_o_M1_level[1, 1:K],
 label = "No migraiton across the North and the South", c=:black,
        line = :dash)
savefig("output/figures/counterfactual/income_no_mig_4060_o_union_conf_bilateral_vs_unilateral_stay.pdf")

plot(years - fill(10.0, K), V_no_conf_to_union_4060_o_M1_level_youngest_ratio[21, 1:K], xlabel = "cohort",
ylabel = "initial expected value", label = "Mississippi", legend = :bottomright, c=:black)
plot!(years - fill(10.0, K), V_no_conf_to_union_4060_o_M1_level_youngest_ratio[10, 1:K], label = "Georgia", c=:black, line = :dash)
plot!(years - fill(10.0, K), V_no_conf_to_union_4060_o_M1_level_youngest_ratio[11, 1:K], label = "Illinois", c=:blue, line = :solid)
plot!(years - fill(10.0, K), V_no_conf_to_union_4060_o_M1_level_youngest_ratio[25, 1:K], label = "New York", c=:blue, line = :dash)
savefig("output/figures/counterfactual/expected_value_cohort_V_no_conf_to_union_4060_o_stay.pdf")

plot(years - fill(10.0, K), V_no_conf_to_union_4060_o_M1_level_youngest_ratio[N+21, 1:K], xlabel = "cohort",
ylabel = "initial expected value", label = "Mississippi: African Americans", legend = :bottomright, c=:black)
plot!(years - fill(10.0, K), V_no_conf_to_union_4060_o_M1_level_youngest_ratio[21, 1:K], label = "Mississippi: Others", c=:black, line = :dash)
plot!(years - fill(10.0, K), V_no_conf_to_union_4060_o_M1_level_youngest_ratio[N+11, 1:K], label = "Illinois: African Americans", c=:blue, line = :solid)
plot!(years - fill(10.0, K), V_no_conf_to_union_4060_o_M1_level_youngest_ratio[11, 1:K], label = "Illinois: Others", c=:blue, line = :dash)
savefig("output/figures/counterfactual/expected_value_cohort_V_no_conf_to_union_4060_o_bo_stay.pdf")

plot(years - fill(10.0, K), V_nounionconf4060_o_M1_level_youngest_ratio[21, 1:K], xlabel = "cohort",
ylabel = "initial expected value", label = "Mississippi", legend = :bottomright, c=:black)
plot!(years - fill(10.0, K), V_nounionconf4060_o_M1_level_youngest_ratio[10, 1:K], label = "Georgia", c=:black, line = :dash)
plot!(years - fill(10.0, K), V_nounionconf4060_o_M1_level_youngest_ratio[11, 1:K], label = "Illinois", c=:blue, line = :solid)
plot!(years - fill(10.0, K), V_nounionconf4060_o_M1_level_youngest_ratio[25, 1:K], label = "New York", c=:blue, line = :dash)
savefig("output/figures/counterfactual/expected_value_cohort_V_nounionconf4060_o_stay.pdf")

plot(years - fill(10.0, K), V_nounionconf4060_o_M1_level_youngest_ratio[N+21, 1:K], xlabel = "cohort",
ylabel = "initial expected value", label = "Mississippi: African Americans", legend = :bottomright, c=:black)
plot!(years - fill(10.0, K), V_nounionconf4060_o_M1_level_youngest_ratio[21, 1:K],
label = "Mississippi: Others", c=:black, line = :dash)
plot!(years - fill(10.0, K), V_nounionconf4060_o_M1_level_youngest_ratio[N+11, 1:K], label = "Illinois: African Americans", c=:blue, line = :solid)
plot!(years - fill(10.0, K), V_nounionconf4060_o_M1_level_youngest_ratio[11, 1:K], label = "Illinois: Others", c=:blue, line = :dash)
savefig("output/figures/counterfactual/expected_value_cohort_V_nounionconf4060_o_bo_stay.pdf")

plot(years - fill(10.0, K), V_no_union_to_conf_4060_o_M1_level_youngest_ratio[21, 1:K], xlabel = "cohort",
ylabel = "initial expected value", label = "Mississippi", legend = :bottomright, c=:black)
plot!(years - fill(10.0, K), V_no_union_to_conf_4060_o_M1_level_youngest_ratio[10, 1:K], label = "Georgia", c=:black, line = :dash)
plot!(years - fill(10.0, K), V_no_union_to_conf_4060_o_M1_level_youngest_ratio[11, 1:K], label = "Illinois", c=:blue, line = :solid)
plot!(years - fill(10.0, K), V_no_union_to_conf_4060_o_M1_level_youngest_ratio[25, 1:K], label = "New York", c=:blue, line = :dash)
savefig("output/figures/counterfactual/expected_value_cohort_V_no_union_to_conf_4060_o_stay.pdf")

plot(years - fill(10.0, K), V_no_union_to_conf_4060_o_M1_level_youngest_ratio[N+21, 1:K], xlabel = "cohort",
ylabel = "initial expected value", label = "Mississippi: African Americans", legend = :bottomright, c=:black)
plot!(years - fill(10.0, K), V_no_union_to_conf_4060_o_M1_level_youngest_ratio[21, 1:K], 
label = "Mississippi: Others", c=:black, line = :dash)
plot!(years - fill(10.0, K), V_no_union_to_conf_4060_o_M1_level_youngest_ratio[N+11, 1:K],
label = "Illinois: African Americans", c=:blue, line = :solid)
plot!(years - fill(10.0, K), V_no_union_to_conf_4060_o_M1_level_youngest_ratio[11, 1:K],
label = "Illinois: Others", c=:blue, line = :dash)
savefig("output/figures/counterfactual/expected_value_cohort_V_no_union_to_conf_4060_o_bo_stay.pdf")

# Confederacy share
plot(years, share_conf_dt[1, 1:K]) 
plot!(years, share_conf_M1_level[1, 1:K])

plot!(years, share_conf_dt[2, 1:K]) 
plot!(years, share_conf_M1_level[2, 1:K])

# increasing migration costs for African Americans from 1940 to 1960 

grids = 21

# store 1970 GDPs 
compare1970gdp = zeros(grids)

C_τ*N*N

path_temp = zeros(size(path1_M1_level))


for i in 1:grids 
        τ_temp = copy(τ)
        τ_temp[(C_τ*N*N+1):(R*C_τ*N*N), 1:3] = 2*(i-1)/(grids-1) * τ[(C_τ*N*N+1):(R*C_τ*N*N), 1:3]

        path_temp = ss_tp(R, C, N, T, λ, λ_2, L_in, s, ν, τ_temp, M1, σ_0_level, σ_1_level, κ_0_level, κ_1_level, A_level, r_bar, η, γ, B, tol, tol_2, maxit, α)

        gdp_temp = GDP(path_temp, T)

        gdp_temp_ratio = gdp_temp ./ gdp1_M1_level

        compare1970gdp[i] = gdp_temp_ratio[1, 4]
end



# store 1970 GDPs 
compare1970gdp_o = zeros(grids)

path_temp = zeros(size(path1_M1_level))


for i in 1:grids 
        τ_temp = copy(τ)
        τ_temp[1:C_τ*N*N, 1:3] =  2*(i-1)/(grids-1) * τ[1:C_τ*N*N, 1:3]

        path_temp = ss_tp(R, C, N, T, λ, λ_2, L_in, s, ν, τ_temp, M1, σ_0_level, σ_1_level, κ_0_level, κ_1_level, A_level, r_bar, η, γ, B, tol, tol_2, maxit, α)

        gdp_temp = GDP(path_temp, T)

        gdp_temp_ratio = gdp_temp ./ gdp1_M1_level

        compare1970gdp_o[i] = gdp_temp_ratio[1, 4]
end

xaxis =  2 * (0:1:(grids-1)) ./ fill((grids-1), grids)

length(xaxis)

plot(xaxis[2:grids], compare1970gdp[2:grids], c = :black, label = "African Americans",
xlab = "migration costs relative to the baseline", ylab = "1970 US income relative to the baseline",
legend = :bottomright)
plot!(xaxis[2:grids], compare1970gdp_o[2:grids], c = :black, line = :dash,
label = "Others")
savefig("output/figures/counterfactual/mig_cost_40_60_income_70.pdf")

# forever 1940
T2 = 101
A_1940 = repeat(A_level[:, 1], 1, T2)
κ_0_1940 = repeat(κ_0_level[:, 1], 1, T2)
κ_1_1940 = repeat(κ_1_level[:, 1], 1, T2)

α_nb_1940 = 1/(s_mat[1, 1] * s_mat[2, 1])
α_b_1940 = 1/(s_mat[8, 1] * s_mat[9, 1])

α_period_1940_p = vcat(repeat([0, 0, α_nb_1940, 0, 0, 0, 0, 0], N), repeat([0, 0, α_b_1940, 0, 0, 0, 0, 0], N))

α_mat_1940 = repeat(α_period_1940_p, 1, T2)

s_1940 = repeat(s[:,1], 1, T2)

B_1940 = repeat(B[:, 1], 1, T2)

τ_1940 = repeat(τ[:, 1], 1, T2)

L_in_1940 = repeat(L_in[:, 1], 1, T2)

path_1940 = ss_tp(R, C, N, T2, λ, λ_2, L_in_1940, s_1940, ν[:,1:T2], τ_1940, M0[:,1:T2], σ_0_level[:,1:T2], 
σ_1_level[:,1:T2], κ_0_1940, κ_1_1940, A_1940, r_bar[:,1:T2], η[:,1:T2], γ[:,1:T2], B_1940, tol, tol_2, maxit, α_mat_1940)

V_1940 = path_1940[:, 3:2+T2]
L_1940 = path_1940[:, (3+T2):(2+2*T2)]
nom_wage_1940 =  path_1940[1:N*R*(C-1), (2+3*T2):(1+4*T2-1)]

gdp_1940 = GDP(path_1940, T2)
local_gdp_1940 = local_GDP(path_1940, N, T2)

share_conf_1940 = share_south(L_1940, union_or_conf, T2)

plot(years[1:5], share_conf_1940[2, 1:5], legend = :topright, label = "African Americans", c = :black,
xlab = "year", ylab = "fraction in the South")
plot!(years[1:5], share_conf_1940[1, 1:5], legend = :topright, label = "Others", c = :black, 
xlab = "year", ylab = "fraction in the South" , line=:dash)
savefig("output/figures/counterfactual/forever1940_conf_share_stay.pdf")

plot(years[1:5], share_conf_1940[2, 1:5], legend = :topright, label = "forever 1940", c = :black,
xlab = "year", ylab = "fraction of African Americans in the South")
plot!(years[1:5], share_conf_dt[2, 1:5], label = "data", c=:blue)
plot!(years[1:5], share_conf_M1_level[2, 1:5], label = "baseline equilibrium", c=:blue, line=:dash)
savefig("output/figures/counterfactual/forever1940_conf_share_vs_baseline_data.pdf")


# forever 1950

A_1950 = repeat(A_level[:, 2], 1, T2)
κ_0_1950 = repeat(κ_0_level[:, 2], 1, T2)
κ_1_1950 = repeat(κ_1_level[:, 2], 1, T2)

α_nb_1950 = 1/(s_mat[1, 2] * s_mat[2, 2])
α_b_1950 = 1/(s_mat[8, 2] * s_mat[9, 2])

α_period_1950_p = vcat(repeat([0, 0, α_nb_1950, 0, 0, 0, 0, 0], N), repeat([0, 0, α_b_1950, 0, 0, 0, 0, 0], N))

α_mat_1950 = repeat(α_period_1950_p, 1, T2)

s_1950 = repeat(s[:,2], 1, T2)

B_1950 = repeat(B[:, 2], 1, T2)

τ_1950 = repeat(τ[:, 2], 1, T2)

L_in_1950 = repeat(L_in[:, 2], 1, T2)

path_1950 = ss_tp(R, C, N, T2, λ, λ_2, L_in_1950, s_1950, ν[:,1:T2], τ_1950, M0[:,1:T2], σ_0_level[:,1:T2], 
σ_1_level[:,1:T2], κ_0_1950, κ_1_1950, A_1950, r_bar[:,1:T2], η[:,1:T2], γ[:,1:T2], B_1950, tol, tol_2, maxit, α_mat_1950)

V_1950 = path_1950[:, 3:2+T2]
L_1950 = path_1950[:, (3+T2):(2+2*T2)]
nom_wage_1950 =  path_1950[1:N*R*(C-1), (2+3*T2):(1+4*T2-1)]

gdp_1950 = GDP(path_1950, T2)
local_gdp_1950 = local_GDP(path_1950, N, T2)

share_conf_1950 = share_south(L_1950, union_or_conf, T2)

# forever 1960
A_1960 = repeat(A_level[:, 3], 1, T2)
κ_0_1960 = repeat(κ_0_level[:, 3], 1, T2)
κ_1_1960 = repeat(κ_1_level[:, 3], 1, T2)

α_nb_1960 = 1/(s_mat[1, 3] * s_mat[2, 3])
α_b_1960 = 1/(s_mat[8, 3] * s_mat[9, 3])

α_period_1960_p = vcat(repeat([0, 0, α_nb_1960, 0, 0, 0, 0, 0], N), repeat([0, 0, α_b_1960, 0, 0, 0, 0, 0], N))

α_mat_1960 = repeat(α_period_1960_p, 1, T2)

s_1960 = repeat(s[:, 3], 1, T2)

B_1960 = repeat(B[:, 3], 1, T2)

τ_1960 = repeat(τ[:, 3], 1, T2)

L_in_1960 = repeat(L_in[:, 3], 1, T2)

path_1960 = ss_tp(R, C, N, T2, λ, λ_2, L_in_1960, s_1960, ν[:,1:T2], τ_1960, M0[:,1:T2], σ_0_level[:,1:T2], 
σ_1_level[:,1:T2], κ_0_1960, κ_1_1960, A_1960, r_bar[:,1:T2], η[:,1:T2], γ[:,1:T2], B_1960, tol, tol_2, maxit, α_mat_1960)

V_1960 = path_1960[:, 3:2+T2]
L_1960 = path_1960[:, (3+T2):(2+2*T2)]
nom_wage_1960 =  path_1960[1:N*R*(C-1), (2+3*T2):(1+4*T2-1)]

gdp_1960 = GDP(path_1960, T2)
local_gdp_1960 = local_GDP(path_1960, N, T2)

share_conf_1960 = share_south(L_1960, union_or_conf, T2)

# Head-Ries style migration costs 

τ_hr_data = CSV.read("../migration_flow/state/output/csv/head_ries_matrix_imp.csv")
τ_hr = convert(Matrix, τ_hr_data)
τ_hr = hcat(τ_hr, repeat(τ_hr[:, K-1], 1, T-(K-1)))

path1_M1_level_hr = ss_tp(R, C, N, T, λ, λ_2, L_in, s, ν, τ_hr, M1, σ_0_level, σ_1_level, κ_0_level, κ_1_level, A_level, r_bar, η, γ, B, tol, tol_2, maxit, α)

gdp1_M1_level_hr = GDP(path1_M1_level_hr, T)

L1_M1_level_hr = path1_M1_level_hr[:, (3+T):(2+2*T)]

local_pop_M1_level_hr = local_pop(L1_M1_level_hr, T, R, N, C)

share_conf_M1_level_hr = share_south(L1_M1_level_hr, union_or_conf, T)



# from 1940 to 1960, African Americans cannot move between the Union and the confederacy
hr_nounionconf4060_b_data = CSV.read("../migration_flow/state/output/csv/head_ries_nounionconf4060_b_matrix.csv")

hr_nounionconf4060_b = make_τ_mat(hr_nounionconf4060_b_data, K_τ-1, C_τ)

path_nounionconf4060_b_M1_level_hr = ss_tp(R, C, N, T, λ, λ_2, L_in, s, ν, hr_nounionconf4060_b, M1, σ_0_level, σ_1_level, κ_0_level, κ_1_level, A_level, r_bar, η, γ, B_d, tol, tol_2, maxit, α)

V_nounionconf4060_b_M1_level_hr = path_nounionconf4060_b_M1_level_hr[:, 3:2+T]

gdp_nounionconf4060_b_M1_level_hr = GDP(path_nounionconf4060_b_M1_level_hr, T)

gdp_ratio_nounionconf4060_b_M1_level_hr = gdp_nounionconf4060_b_M1_level_hr ./ gdp1_M1_level_hr_B_d



# compare US labor incomes 
plot(years, US_income_dt[1, 1:K], legend = :topleft, label = "data", c = :blue,
xlab = "year", ylab = "US labor income")
plot!(years, gdp1_M1_level_hr[1, 1:K], label = "Head-Ries", c=:black, line=:dash)
plot!(years, gdp1_M1_level[1, 1:K], label = "fixed effects", c=:black)
# savefig("output/figures/model_fit/income_model_data_hr_B_d_hikaku.pdf")

# compare population vectors at location-age-race levels
pop_cor_hr = zeros(K)
for i in 1:K
        pop_cor_hr[i] = cor(L_in[:, i], L1_M1_level_hr[:, i])
end

plot(years, pop_cor, label = "fixed effects", xlab = "year", ylab = "correlation",
c = :black, ylims = [0.95, 1.02])
plot!(years, pop_cor_hr, label = "Head-Ries", xlab = "year", ylab = "correlation",
c = :black, line=:dash, ylims = [0.95, 1.02])
# savefig("output/figures/model_fit/correlation_pop_model_data_hr_B_d_hikaku.pdf")

# compare populations vectors at location-race lavels 
pop_cor_loc_hr_B_d = zeros(K)
for i in 1:K
        pop_cor_loc_hr_B_d[i] = cor(local_pop_dt[:, i], local_pop_M1_level_hr_B_d[:, i])
end

plot(years, pop_cor_loc, label = "fixed effects", xlab = "year", ylab = "correlation",
c = :black, ylims = [0.85, 1.02])
plot!(years, pop_cor_loc_hr_B_d, label = "Head-Ries", xlab = "year", ylab = "correlation",
c = :black, line=:dash, ylims = [0.85, 1.02])
# savefig("output/figures/model_fit/correlation_pop_model_data_loc_hikaku.pdf")

plot(years, share_conf_dt[1, 1:K]) 
plot!(years, share_conf_M1_level[1, 1:K])
plot!(years, share_conf_M1_level_hr[1, 1:K])

plot!(years, share_conf_dt[2, 1:K]) 
plot!(years, share_conf_M1_level[2, 1:K])
plot!(years, share_conf_M1_level_hr[2, 1:K])

# IV 2

# elasticity of substitution between ages 
σ_0_level_iv2 = fill(1/0.5616, 1, T) 

# elasticity of substitution between races (within ages)
σ_1_level_iv2 = fill(1/0.1246, 1, T) 

# productivity
A_data_level_iv2 = CSV.read("../elas_sub/output/csv/loc_prod_matrix_level_iv2.csv")

A_mat_1_level_iv2 =  convert(Matrix, A_data_level_iv2)

A_mat_2_level_iv2 = repeat(A_mat_1_level_iv2[:, K], 1, (T-K))

A_level_iv2 = hcat(A_mat_1_level_iv2, A_mat_2_level_iv2)

# cohort-specific productivity 
κ_0_data_level_iv2 = CSV.read("../elas_sub/output/csv/age_prod_matrix_level_iv2.csv")
κ_0_mat_1_level_iv2 = convert(Matrix, κ_0_data_level_iv2)

κ_0_mat_2_level_iv2 = repeat(κ_0_mat_1_level_iv2[:, K], 1, (T-K))

κ_0_level_iv2 = hcat(κ_0_mat_1_level_iv2, κ_0_mat_2_level_iv2)

# race-specific productivity (within cohorts)
κ_1_data_level_iv2 = CSV.read("../elas_sub/output/csv/race_prod_matrix_level_iv2.csv")
κ_1_mat_1_level_iv2 = convert(Matrix, κ_1_data_level_iv2)

κ_1_mat_2_level_iv2 = repeat(κ_1_mat_1_level_iv2[:, K], 1, (T-K))

κ_1_level_iv2 = hcat(κ_1_mat_1_level_iv2, κ_1_mat_2_level_iv2)

A_level_iv2
A_level

κ_0_level_iv2
κ_0_level

κ_1_level_iv2
κ_1_level

# the baseline equilibrium
path1_M1_level_iv2 = ss_tp(R, C, N, T, λ, λ_2, L_in, s, ν, τ, M1, σ_0_level_iv2, σ_1_level_iv2,  κ_0_level_iv2, κ_1_level_iv2, A_level_iv2, r_bar, η, γ, B, tol, tol_2, maxit, α)
save_path1_M1_level_iv2 = convert(DataFrame, path1_M1_level_iv2)
CSV.write("output/csv/path1_M1_level_iv2.csv", save_path1_M1_level_iv2)




V1_M1_level_iv2 = path1_M1_level_iv2[:, 3:2+T]
L1_M1_level_iv2 = path1_M1_level_iv2[:, (3+T):(2+2*T)]
nom_wage1_M1_level_iv2 =  path1_M1_level_iv2[1:N*R*(C-1), (2+3*T):(1+4*T-1)]

gdp1_M1_level_iv2 = GDP(path1_M1_level_iv2, T)
local_gdp1_M1_level_iv2 = local_GDP(path1_M1_level_iv2, N, T)

share_south_M1_level_iv2 = share_south(L1_M1_level_iv2, north_or_south, T)
share_conf_M1_level_iv2 = share_south(L1_M1_level_iv2, union_or_conf, T)

local_pop_M1_level_iv2 = local_pop(L1_M1_level_iv2, T, R, N, C)

V_youngest1_M1_level_iv2 = get_V_youngest(V1_M1_level_iv2)

# African Americans cannot move across the Union and the Confederacy from 1940 to 1960
path_nounionconf_b_M1_level_iv2 = ss_tp(R, C, N, T, λ, λ_2, L_in, s, ν, τ_nounionconf4060_b, M1, σ_0_level_iv2, σ_1_level_iv2,
 κ_0_level_iv2, κ_1_level_iv2, A_level_iv2, r_bar, η, γ, B, tol, tol_2, maxit, α)

save_path_nounionconf_b_M1_level_iv2 = convert(DataFrame, path_nounionconf_b_M1_level_iv2)
CSV.write("output/csv/path_nounionconf_b_M1_level_iv2.csv", save_path_nounionconf_b_M1_level_iv2)

V_nounionconf_b_M1_level_iv2 = path_nounionconf_b_M1_level_iv2[:, 3:2+T]

gdp_nounionconf_b_M1_level_iv2 = GDP(path_nounionconf_b_M1_level_iv2, T)

gdp_ratio_nounionconf_b_M1_level_iv2 = gdp_nounionconf_b_M1_level_iv2 ./ gdp1_M1_level_iv2

gdp1_M1_level_iv2 ./ gdp_nounionconf_b_M1_level_iv2

V_nounionconf_b_M1_level_iv2_youngest = get_V_youngest(V_nounionconf_b_M1_level_iv2)
V_nounionconf_b_M1_level_iv2_youngest_ratio = V_nounionconf_b_M1_level_iv2_youngest ./ V_youngest1_M1_level_iv2

save_V_no_union_to_conf_4060_o_M1_level_youngest_ratio = hcat(V_no_union_to_conf_4060_o_M1_level_youngest_ratio[1:N, 1],
V_no_union_to_conf_4060_o_M1_level_youngest_ratio[(N+1):R*N, 1])