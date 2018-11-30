############################################
### Neural Network with One Hidden Layer ###
############################################

## Packages
using Distributions
using Plots
pyplot()

########## Generate Some Data
## Define Circles
cx1 = 0.3
cy1 = 0.3
r1 = 0.3

cx2 = -0.5
cy2 = -0.5
r2 = 0.3

## Generate data points
nₛ = 400
X_data = 2*(rand(nₛ,2) .-0.5)

## Check if Data is Inside Circle
C1 = (X_data[:,1] .-cx1).^2 .+ (X_data[:,2] .-cy1).^2 .<r1^2
C2 = (X_data[:,1] .-cx2).^2 .+ (X_data[:,2] .-cy2).^2 .<r2^2
T = 1.0*(C1.|C2)

## Plotting
scatter(X_data[:,1], X_data[:,2], legend=false)
scatter(X_data[:,1], X_data[:,2], zcolor=T, colorbar=false, legend=false)

########## X_tilde
X_tilde = hcat(X_data,ones(nₛ))
nₓ = 2
nₕ = 6

########## Define Functions
lfun(z)=1 ./(1 .+exp.(-z))

## Prediction Function
function Prediction(Xtil, w, v)
    H = lfun(Xtil*w)
    H_tilde = hcat(H, ones(length(H[:,1])))
    P = lfun(H_tilde*v)

    return P
end

## Gradient Functions
function derVW(Xtil, P, T, w, v)
    H = lfun(Xtil*w)
    nₛ = length(H[:,1])
    H_tilde = hcat(H, ones(nₛ))

    Δᵖ = P-T
    dv = (H_tilde'*Δᵖ)./nₛ

    Θ = Δᵖ*v[1:end-1]'
    Δʰ = (H.*(1 .-H)).*Θ
    dw = (Xtil'*Δʰ)./nₛ
    return dv, dw
end

## Cost Function
function Cost(X, P, T)
    n = length(T)
    return sum(T.*log.(P) .+(1 .-T).*log.(1 .-P)) /(-n)
end

##### Testing #####
P = Prediction(X_tilde, w, v)
dv, dw = derVW(X_tilde, P, T, w, v)

########## Optimization Loop
function OptimizationLoop(iter, α, X_tilde, T, w_vec, v_vec, Cost_vec)
    for i in 2:iter
        P = Prediction(X_tilde, w_vec, v_vec)
        v_up, w_up = derVW(X_tilde, P, T, w_vec, v_vec)
        w_vec -= α.*w_up
        v_vec -= α.*v_up
        Cost_vec[i] = Cost(P, P, T)
    end

    return P
end

iter = 1000
α = 0.1
w_vec = 0.1.*ones(nₓ+1,nₕ)
v_vec = 0.1.*ones(nₕ+1,1)
Cost_vec = zeros(iter)
P₀ = Prediction(X_tilde, w_vec, v_vec)
Cost_vec[1] = Cost(X_tilde, P₀, T)

Pred = OptimizationLoop(iter, α, X_tilde, T, w_vec, v_vec, Cost_vec)

########## Plotting
plot(collect(1:1:iter), Cost_vec, xscale=:log10, title="Cost", xlabel="iterations", ylabel="Cost", legend=false)

for i in 1:length(Pred)
    Pred[i] < 0.5 ? Pred[i] = 0 : Pred[i] = 1
end
scatter(X_data[:,1], X_data[:,2], zcolor=Pred, colorbar=false, legend=false)
