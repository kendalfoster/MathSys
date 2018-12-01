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
σϵ = 2
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
plot_52i = plot(X, title="Autoregressive Model Simulation", xlabel="Timestep", ylabel="X Value", legend=false)
savefig(plot_52i, "Images/plot_52i.png")

## Calculate the Sample Mean and Sample Variance
X_mean = sum(X)/N
X_var = sum((X.-X_mean).^2)/N

## Theoretical Mean and Variance
μ_th = c/(1-ϕ)
σ_th = (σϵ^2)/(1-ϕ^2)

## Percent Error
μ_pe = 100*abs(μ_th-X_mean)/μ_th
σ_pe = 100*abs(σ_th-X_var)/σ_th

## Theoretical Autocovariance
L = 10 # [-L,L] is range for autocorrelation
ac_th = zeros(2*L+1)
for n in 1:2*L+1
    ac_th[n] = (σϵ^2/(1-ϕ^2))*ϕ^abs(n-L-1)
end

## Empirical Autocovariance
covarp = zeros(N)
covarp[1] = var(X)
covarn = zeros(N)
X_r = reverse(X)

for j =2:N
    K = (N-j-1)
    for k=1:K
        n = j-1
        k_= k+n
        covarp[j] += ((X[k]-X_mean)*(X[k_]-X_mean))
        covarn[j] += ((X_r[k]-X_mean)*(X_r[k_]-X_mean))
    end
    covarp[j] = covarp[j]/K
    covarn[j] = covarn[j]/K
end

covarn = reverse(covarn)
covar = vcat(covarn[1:end-1],covarp)

## Plotting
plot_52ii = plot(collect(-L:1:L), ac_th, title="Theoretical and Empirical Autocorrelation", xlabel="Value of n", ylabel="Autocorrelation", label="Theoretical")
plot!(collect(-L:1:L), covar[N-L:N+L], label="Empirical")
savefig(plot_52ii, "Images/plot_52ii")
