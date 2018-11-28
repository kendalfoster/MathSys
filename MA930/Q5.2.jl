#################################################
################ Question 5.2 ###################
#################################################

## Packages
using Distributions, Random

## Use Julia Plots, with PyPlot backend
using Plots
pyplot()

## Define Parameters
c = 1
ϕ = 0.5
σ_ϵ = 2
N = 5000 # timesteps
st = sqrt(2)/10

## Make vectors
X = st*ones(N)
ϵ = rand(Normal(0,σϵ),N)

## Define the model
for n in 2:N
    X[n] = c + ϕ*X[n-1] + ϵ[n]
end

## Plotting
plot(X, title="Autoregressive Model Simulation", xlabel="Timestep", ylabel="X Value", legend=false)

## Calculate the Sample Mean and Sample Variance
X_mean = sum(X)/N
X_var = sum((X.-X_mean).^2)/N

## Theoretical Mean and Variance
μ_th = c/(1-ϕ)
σ_th = (σϵ^2)/(1-ϕ^2)

## Theoretical Autocorrelation
L = 10 # [-L,L] is range for autocorrelation
ac_th = zeros(2*L+1)
for n in 1:2*L+1
    ac_th[n] = (σ_ϵ^2/(1-ϕ^2))*ϕ^abs(n-L-1)
end

plot(collect(-L:1:L), ac_th, title="Theoretical and Computational Autocorrelation", xlabel="Value of n", ylabel="Autocorrelation", label="Theoretical")
