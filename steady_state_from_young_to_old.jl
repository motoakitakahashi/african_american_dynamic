# Motoaki Takahashi

# Julia 1.7.2

# February 2022

# I compute a steady state of the model

# for the main desktop
cd("D:/onedrive/OneDrive - The Pennsylvania State University/dynamic/code")

# for the laptop
# cd("C:/Users/takah/OneDrive - The Pennsylvania State University/dynamic/code")

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

# the number of places 
N = 2

# the number of cohorts in a given period
C = 8

# the number of races 
R = 2

# migration elesticity
ν = 2.0 # this follows Caliendo, Opromolla, Parro, and Sforza

# elasticity of substitution between ages 
σ_0 = 5.0 # this follows Ottaviano and Peri, Manacorda et al, and Yuta's JMP

# elasticity of substitution between races (within ages)
σ_1 = 11.0 # this follows Boustan

# Cobb-Douglas share of housing
γ = 0.25


# amenities 
B = [   1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
# the amenity for the youngest race 1 in place 1, ..., the amenity for the oldest race 1 in place 1, ...
# the amenity for the youngest race 1 in place 2, ..., the amenity for the oldest race 1 in place 2, ...
# the amenity for the youngest race 2 in place 1, ..., the amenity for the oldest race 2 in place 1, ...
# the amenity for the youngest race 2 in place 2, ..., the amenity for the oldest race 2 in place 2, ...

# survival probabilities
s = [   1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
# the survival probability for the youngest race 1, ..., the survival probability for the second oldest race 1, ...
# the survival probability for the youngest race 2, ..., the survival probability for the second oldest race 2, ...


# migration costs 
# note that the oldest can no longer migrate

τ = [   0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
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

# place-specific shifter of rent 
r_bar = [0.1, 0.1]

# place-specific elasticity of rent 
η = [0.5, 0.5]





# productivity
A = [1.0, 1.0]
# the productivity in place 1, the productivity in place 2

# cohort-specific productivity 
κ_0 = [1/7, 1/7, 1/7, 1/7, 1/7, 1/7, 1/7,
        1/7, 1/7, 1/7, 1/7, 1/7, 1/7, 1/7]
# the relative productivity of the second youngest in place 1, ..., the relative productivity of the oldest in place 1,
# the relative productivity of the second youngest in place 2, ..., the relative productivity of the oldest in place 2

# race-specific productivity (within cohorts)
κ_1 = [1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2,
        1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2,
        1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2,
        1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2]
# the relative productivity of race1 within the second youngest in place 1, ..., the relative productivity of race1 within the oldest in place 1,
# the relative productivity of race1 within the second youngest in place 2, ..., the relative productivity of race1 within the oldest in place 2,
# the relative productivity of race2 within the second youngest in place 1, ..., the relative productivity of race2 within the oldest in place 1,
# the relative productivity of race2 within the second youngest in place 2, ..., the relative productivity of race2 within the oldest in place 2

# fertility per cohort-race-place 
# I guess I need to restrict the value of α to get a steady state, 
α = [   0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,]
# the fertility of the youngest race 1 in place 1, ..., the fertility of the oldest race 1 in place 1,
# the fertility of the youngest race 2 in place 1, ..., the fertility of the oldest race 2 in place 1,
# the fertility of the youngest race 2 in place 1, ..., the fertility of the oldest race 2 in place 1,
# the fertility of the youngest race 2 in place 2, ..., the fertility of the oldest race 2 in place 2


# function that maps wages, rents, and amenties to period utility

# rents 
r = [1.0, 1.0] # rent in place 1, rent in place 2

# wages 
w = [   10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0,
        10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0,
        10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0,
        10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0]
# wage of the second youngest race 1 in place 1, ..., wage of the oldest race 1 in place 1,
# wage of the second youngest race 1 in place 2, ..., wage of the oldest race 1 in place 2,
# wage of the second youngest race 2 in place 1, ..., wage of the oldest race 2 in place 1,
# wage of the second youngest race 2 in place 2, ..., wage of the oldest race 2 in place 2



function utility(w, r, B, R, N, C, γ)
        temp = log.(w./ kron(r .^ γ, fill(1.0, ((C-1)*R))))
        u = []

        for i in 1:(R*N)
                u = [u; 0.0; temp[((i-1)*(C-1)+1):(i*(C-1))] ]
        end

        u = u + log.(B)

        return u
end

u = utility(w, r, B, R, N, C, γ)


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

V = expected_value(u, s, ν, τ, R, N, C)

V_mat = reshape(V, C, N*R)'


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

μ = migration_rate(s, V, ν, τ, R, N, C)

# for now, let's assign some values to populations
L = [   10.0, 10.0*0.9, 10.0*(0.9)^2, 10.0*(0.9)^3, 10.0*(0.9)^4, 10.0*(0.9)^5, 10.0*(0.9)^6, 10.0*(0.9)^7,
        10.0, 10.0*0.9, 10.0*(0.9)^2, 10.0*(0.9)^3, 10.0*(0.9)^4, 10.0*(0.9)^5, 10.0*(0.9)^6, 10.0*(0.9)^7, 
        10.0, 10.0*0.9, 10.0*(0.9)^2, 10.0*(0.9)^3, 10.0*(0.9)^4, 10.0*(0.9)^5, 10.0*(0.9)^6, 10.0*(0.9)^7,
        10.0, 10.0*0.9, 10.0*(0.9)^2, 10.0*(0.9)^3, 10.0*(0.9)^4, 10.0*(0.9)^5, 10.0*(0.9)^6, 10.0*(0.9)^7]
# population of the youngest race 1 in place 1, ..., population of the oldest race 1 in place 1,
# population of the youngest race 1 in place 2, ..., population of the oldest race 1 in place 2,
# population of the youngest race 2 in place 1, ..., population of the oldest race 2 in place 1,
# population of the youngest race 2 in place 2, ..., population of the oldest race 2 in place 2.


# here, for computation of steady states, we ignore immigration from abroad
# the youngest cohort is pinned down by fertility, the other cohorts are pinned down by migration flows (and death/out-migration to abroad)
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

L_new = population(μ, s, L, α, C, R, N)

# the following function aggregates race-cohort level labor to cohort-level labor
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

L_cohort = agg_labor_cohort(L, σ_1, κ_1, N, R, C)

# the following function aggregates race-cohort level labor to cohort-level labor
function agg_labor_place(L_cohort, σ_0, κ_0, N, R, C)
        L_cohort_mat = reshape(L_cohort, (C-1), N)'
        κ_0_mat = reshape(κ_0, (C-1), N)'
        
        temp = κ_0_mat .^ (1/σ_0) .* L_cohort_mat .^ ((σ_0-1)/σ_0)
        L_place_mat = sum(temp, dims = 2) .^ (σ_0/(σ_0-1))

        return L_place_mat

end

L_place = agg_labor_place(L_cohort, σ_0, κ_0, N, R, C)





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


w = wage(σ_0, σ_1, A, κ_0, κ_1, N, R, C, L_place, L_cohort, L)



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

r = rent(w, L, r_bar, η, γ, N, R, C)


tol = 10 ^ (-7)
maxit = 1000

# the dumpening parameter for iteration
λ = 0.5


function steady_state(tol, maxit, N, C, R, ν, σ_0, σ_1, γ, B, s, τ, r_bar, η, A, κ_0, κ_1, α, L_in, λ)
        L = L_in 

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

        output = [L fill(dif, size(L)) fill(count, size(L))]

        return output
end

L_in = [10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0,
        10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 
        10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0,
        10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0]

A = [1.0, 1.0]

L_dif_SS = steady_state(tol, maxit, N, C, R, ν, σ_0, σ_1, γ, B, s, τ, r_bar, η, A, κ_0, κ_1, α, L_in, λ)



# comparative statics (or more precisely comparative steady states)

# comparative statics with respect to productivity

A_1_range = 2.0:0.2:8.0

K = size(A_1_range)[1]

L_SS_cs_A_1 = zeros((R*N*C), K)
w_SS_cs_A_1 = zeros((R*N*(C-1)), K)
r_SS_cs_A_1 = zeros(N, K)
V_SS_cs_A_1 = zeros((R*N*C), K)
μ_SS_cs_A_1 = zeros((R*N*N*(C-1)), K)


dif_SS_cs_A_1 = zeros((R*N*C), K)

for i in 1:K
        A = [A_1_range[i], 5.0]
        temp = steady_state(tol, maxit, N, C, R, ν, σ_0, σ_1, γ, B, s, τ, r_bar, η, A, κ_0, κ_1, α, L_in, λ)
        L = temp[:, 1]
        L_cohort = agg_labor_cohort(L, σ_1, κ_1, N, R, C)
        L_place = agg_labor_place(L_cohort, σ_0, κ_0, N, R, C)

        w = wage(σ_0, σ_1, A, κ_0, κ_1, N, R, C, L_place, L_cohort, L)
        r = rent(w, L, r_bar, η, γ, N, R, C)
        u = utility(w, r, B, R, N, C, γ)
        V = expected_value(u, s, ν, τ, R, N, C)
        μ = migration_rate(s, V, ν, τ, R, N, C)

        L_SS_cs_A_1[:, i] = temp[:, 1]
        dif_SS_cs_A_1[:, i] = temp[:, 2]
        w_SS_cs_A_1[:, i] = w
        r_SS_cs_A_1[:, i] = r
        V_SS_cs_A_1[:, i] = V 
        μ_SS_cs_A_1[:, i] = μ




end

L_SS_cs_A_1_0_1 = reshape(L_SS_cs_A_1[:, 1], C, (R*N))'

L_SS_cs_A_1_1_9 = reshape(L_SS_cs_A_1[:, K], C, (R*N))'

plot(A_1_range, L_SS_cs_A_1[1, :])

plot!(A_1_range, L_SS_cs_A_1[((C+1)), :])

plot!(A_1_range, L_SS_cs_A_1[C, :])

plot!(A_1_range, L_SS_cs_A_1[2*C, :])

plot(A_1_range, L_SS_cs_A_1[2*C+1, :])
plot!(A_1_range, L_SS_cs_A_1[2*C+C, :])
plot!(A_1_range, L_SS_cs_A_1[3*C+1, :])
plot!(A_1_range, L_SS_cs_A_1[4*C, :])

plot(A_1_range, r_SS_cs_A_1[1, :])
plot!(A_1_range, r_SS_cs_A_1[2, :])

plot(A_1_range, w_SS_cs_A_1[1, :])
plot!(A_1_range, w_SS_cs_A_1[(C-1)+1, :])

plot(A_1_range, w_SS_cs_A_1[(C-1), :])
plot!(A_1_range, w_SS_cs_A_1[2*(C-1), :])



plot(A_1_range, V_SS_cs_A_1[1, :])
plot!(A_1_range, V_SS_cs_A_1[2, :])
plot!(A_1_range, V_SS_cs_A_1[3, :])
plot!(A_1_range, V_SS_cs_A_1[4, :])
plot!(A_1_range, V_SS_cs_A_1[5, :])
plot!(A_1_range, V_SS_cs_A_1[6, :])
plot!(A_1_range, V_SS_cs_A_1[7, :])
plot!(A_1_range, V_SS_cs_A_1[C, :])
plot!(A_1_range, V_SS_cs_A_1[(C+1), :])
plot!(A_1_range, V_SS_cs_A_1[(C+2), :])
plot!(A_1_range, V_SS_cs_A_1[(C+3), :])
plot!(A_1_range, V_SS_cs_A_1[(C+4), :])
plot!(A_1_range, V_SS_cs_A_1[(C+5), :])
plot!(A_1_range, V_SS_cs_A_1[(C+6), :])
plot!(A_1_range, V_SS_cs_A_1[(C+7), :])
plot!(A_1_range, V_SS_cs_A_1[(2*C), :])




# comparative statics with respect to amenity

A = [5.0, 5.0]

B_1_1_range = 0.5:0.05:1.5

K = size(B_1_1_range)[1]

L_SS_cs_B_1_1 = zeros((R*N*C), K)
w_SS_cs_B_1_1 = zeros((R*N*(C-1)), K)
r_SS_cs_B_1_1 = zeros(N, K)
V_SS_cs_B_1_1 = zeros((R*N*C), K)
μ_SS_cs_B_1_1 = zeros((R*N*N*(C-1)), K)


dif_SS_cs_B_1_1 = zeros((R*N*C), K)

for i in 1:K
        B = [fill(B_1_1_range[i], C); fill(1.0, (C*(N+1)))]
        temp = steady_state(tol, maxit, N, C, R, ν, σ_0, σ_1, γ, B, s, τ, r_bar, η, A, κ_0, κ_1, α, L_in, λ)
        L = temp[:, 1]
        L_cohort = agg_labor_cohort(L, σ_1, κ_1, N, R, C)
        L_place = agg_labor_place(L_cohort, σ_0, κ_0, N, R, C)

        w = wage(σ_0, σ_1, A, κ_0, κ_1, N, R, C, L_place, L_cohort, L)
        r = rent(w, L, r_bar, η, γ, N, R, C)
        u = utility(w, r, B, R, N, C, γ)
        V = expected_value(u, s, ν, τ, R, N, C)
        μ = migration_rate(s, V, ν, τ, R, N, C)

        L_SS_cs_B_1_1[:, i] = temp[:, 1]
        dif_SS_cs_B_1_1[:, i] = temp[:, 2]
        w_SS_cs_B_1_1[:, i] = w
        r_SS_cs_B_1_1[:, i] = r
        V_SS_cs_B_1_1[:, i] = V 
        μ_SS_cs_B_1_1[:, i] = μ




end

plot(B_1_1_range, L_SS_cs_B_1_1[1, :])

plot!(B_1_1_range, L_SS_cs_B_1_1[((C+1)), :])

plot!(B_1_1_range, L_SS_cs_B_1_1[C, :])

plot!(B_1_1_range, L_SS_cs_B_1_1[2*C, :])

plot(B_1_1_range, L_SS_cs_B_1_1[2*C+1, :])
plot!(B_1_1_range, L_SS_cs_B_1_1[2*C+C, :])
plot!(B_1_1_range, L_SS_cs_B_1_1[3*C+1, :])
plot!(B_1_1_range, L_SS_cs_B_1_1[4*C, :])

plot(B_1_1_range, r_SS_cs_B_1_1[1, :])
plot!(B_1_1_range, r_SS_cs_B_1_1[2, :])

plot(B_1_1_range, w_SS_cs_B_1_1[1, :])
plot!(B_1_1_range, w_SS_cs_B_1_1[(C-1)+1, :])

plot(B_1_1_range, w_SS_cs_B_1_1[(C-1), :])
plot!(B_1_1_range, w_SS_cs_B_1_1[2*(C-1), :])



plot(B_1_1_range, V_SS_cs_B_1_1[1, :])
plot!(B_1_1_range, V_SS_cs_B_1_1[2, :])
plot!(B_1_1_range, V_SS_cs_B_1_1[3, :])
plot!(B_1_1_range, V_SS_cs_B_1_1[4, :])
plot!(B_1_1_range, V_SS_cs_B_1_1[5, :])
plot!(B_1_1_range, V_SS_cs_B_1_1[6, :])
plot!(B_1_1_range, V_SS_cs_B_1_1[7, :])
plot!(B_1_1_range, V_SS_cs_B_1_1[C, :])
plot!(B_1_1_range, V_SS_cs_B_1_1[(C+1), :])
plot!(B_1_1_range, V_SS_cs_B_1_1[(C+2), :])
plot!(B_1_1_range, V_SS_cs_B_1_1[(C+3), :])
plot!(B_1_1_range, V_SS_cs_B_1_1[(C+4), :])
plot!(B_1_1_range, V_SS_cs_B_1_1[(C+5), :])
plot!(B_1_1_range, V_SS_cs_B_1_1[(C+6), :])
plot!(B_1_1_range, V_SS_cs_B_1_1[(C+7), :])
plot!(B_1_1_range, V_SS_cs_B_1_1[(2*C), :])


plot(B_1_1_range, V_SS_cs_B_1_1[2*C+1, :])
plot!(B_1_1_range, V_SS_cs_B_1_1[2*C+2, :])
plot!(B_1_1_range, V_SS_cs_B_1_1[2*C+3, :])
plot!(B_1_1_range, V_SS_cs_B_1_1[2*C+4, :])
plot!(B_1_1_range, V_SS_cs_B_1_1[2*C+5, :])
plot!(B_1_1_range, V_SS_cs_B_1_1[2*C+6, :])
plot!(B_1_1_range, V_SS_cs_B_1_1[2*C+7, :])
plot(B_1_1_range, V_SS_cs_B_1_1[2*C+C, :])
plot!(B_1_1_range, V_SS_cs_B_1_1[3*C+1, :])
plot!(B_1_1_range, V_SS_cs_B_1_1[3*C+2, :])
plot!(B_1_1_range, V_SS_cs_B_1_1[3*C+3, :])
plot!(B_1_1_range, V_SS_cs_B_1_1[3*C+4, :])
plot!(B_1_1_range, V_SS_cs_B_1_1[3*C+5, :])
plot!(B_1_1_range, V_SS_cs_B_1_1[3*C+6, :])
plot!(B_1_1_range, V_SS_cs_B_1_1[3*C+7, :])
plot!(B_1_1_range, V_SS_cs_B_1_1[3*C+C, :])


# comparative statics with respect to age-specific productivity

A = [5.0, 5.0]

B = [   1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

κ_0_1_range = 1/700:1/700:200/700

K = size(κ_0_1_range)[1]

L_SS_cs_κ_0_1 = zeros((R*N*C), K)
w_SS_cs_κ_0_1 = zeros((R*N*(C-1)), K)
r_SS_cs_κ_0_1 = zeros(N, K)
V_SS_cs_κ_0_1 = zeros((R*N*C), K)
μ_SS_cs_κ_0_1 = zeros((R*N*N*(C-1)), K)


dif_SS_cs_κ_0_1 = zeros((R*N*C), K)

# let's move 30s productivity
for i in 1:K 
        κ_0_1_3 = κ_0_1_range[i]
        κ_0 = [ (1-κ_0_1_3)/7, (1-κ_0_1_3)/7, κ_0_1_3, (1-κ_0_1_3)/7, (1-κ_0_1_3)/7, (1-κ_0_1_3)/7, (1-κ_0_1_3)/7,
                1/7, 1/7, 1/7, 1/7, 1/7, 1/7, 1/7]

        temp = steady_state(tol, maxit, N, C, R, ν, σ_0, σ_1, γ, B, s, τ, r_bar, η, A, κ_0, κ_1, α, L_in, λ)

        L = temp[:, 1]
        L_cohort = agg_labor_cohort(L, σ_1, κ_1, N, R, C)
        L_place = agg_labor_place(L_cohort, σ_0, κ_0, N, R, C)

        w = wage(σ_0, σ_1, A, κ_0, κ_1, N, R, C, L_place, L_cohort, L)
        r = rent(w, L, r_bar, η, γ, N, R, C)
        u = utility(w, r, B, R, N, C, γ)
        V = expected_value(u, s, ν, τ, R, N, C)
        μ = migration_rate(s, V, ν, τ, R, N, C)

        L_SS_cs_κ_0_1[:, i] = L
        dif_SS_cs_κ_0_1[:, i] = temp[:, 2]
        w_SS_cs_κ_0_1[:, i] = w
        r_SS_cs_κ_0_1[:, i] = r
        V_SS_cs_κ_0_1[:, i] = V 
        μ_SS_cs_κ_0_1[:, i] = μ 

        

end

plot(κ_0_1_range, L_SS_cs_κ_0_1[1, :])
plot!(κ_0_1_range, L_SS_cs_κ_0_1[2, :])
plot!(κ_0_1_range, L_SS_cs_κ_0_1[3, :])
plot!(κ_0_1_range, L_SS_cs_κ_0_1[4, :])
plot!(κ_0_1_range, L_SS_cs_κ_0_1[5, :])
plot!(κ_0_1_range, L_SS_cs_κ_0_1[6, :])
plot!(κ_0_1_range, L_SS_cs_κ_0_1[7, :])
plot!(κ_0_1_range, L_SS_cs_κ_0_1[C, :])

plot(κ_0_1_range, L_SS_cs_κ_0_1[(C+1), :])
plot!(κ_0_1_range, L_SS_cs_κ_0_1[(C+2), :])
plot!(κ_0_1_range, L_SS_cs_κ_0_1[(C+3), :])
plot!(κ_0_1_range, L_SS_cs_κ_0_1[(C+4), :])
plot!(κ_0_1_range, L_SS_cs_κ_0_1[(C+5), :])
plot!(κ_0_1_range, L_SS_cs_κ_0_1[(C+6), :])
plot!(κ_0_1_range, L_SS_cs_κ_0_1[(C+7), :])
plot!(κ_0_1_range, L_SS_cs_κ_0_1[2*C, :])

plot(κ_0_1_range, L_SS_cs_κ_0_1[(2*C+1), :])
plot!(κ_0_1_range, L_SS_cs_κ_0_1[(2*C+2), :])
plot!(κ_0_1_range, L_SS_cs_κ_0_1[(2*C+3), :])
plot!(κ_0_1_range, L_SS_cs_κ_0_1[(2*C+4), :])
plot!(κ_0_1_range, L_SS_cs_κ_0_1[(2*C+5), :])
plot!(κ_0_1_range, L_SS_cs_κ_0_1[(2*C+6), :])
plot!(κ_0_1_range, L_SS_cs_κ_0_1[(2*C+7), :])
plot!(κ_0_1_range, L_SS_cs_κ_0_1[3*C, :])


