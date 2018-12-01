#################################################
################ Question 5.1 ###################
#################################################

## Packages
using DelimitedFiles
using DataFrames
using NaNMath
using Statistics

## Use Julia Plots, with PyPlot backend
using Plots
pyplot()

###### Read and clean data
data = readdlm("exoplanet-data.txt")
df = convert(DataFrame, data)
N = length(df[:,2])
for i in 1:N
    if isnan(df[i,2])==true
        df[i,2] = NaNMath.mean([df[i-2,2], df[i-1,2], df[i+1,2], df[i+2,2]])
    end
end
# # Boxcar smoothing
# for i in 1:N
#     if i<=2
#         df[i,2] = mean([3*df[i,2], df[i+1,2], df[i+2,2]])
#     elseif i>=N-2
#         df[i,2] = mean([df[i-2,2], df[i-1,2], 3*df[i,2]])
#     else
#         df[i,2] = mean([df[i-2,2], df[i-1,2], df[i,2], df[i+1,2], df[i+2,2]])
#     end
# end

###### Initial look at data
plot(df[:,1],df[:,2], title="Light Intensity of a Star", xlabel="Time, Earth Days", ylabel="Intensity",legend=false)

###### Part a: Estimate the orbital period of the planet
## Crop the data to determine regular light intensity
# Have to make this a function because Julia is absurd
function reg_intensity(df)
    ε = 5
    dum_df = copy(df)
    for i in 1:5 # number of dips (clear from graph)
        m = argmin(dum_df[:,2])
        println("m = $m")
        dum_df = dum_df[setdiff(1:end, m-ε:m+ε), :] # remove the rows of data around the minimum
    end
    return dum_df
end

function idx_mins(df)
    another = copy(df)
    I₀_vec = reg_intensity(another)
    I₀ = mean(I₀_vec[:,2])

    idx_mins = []
    ε = 5
    dum_df = copy(df)
    for i in 1:5 # number of dips (clear from graph)
        m = argmin(dum_df[:,2])
        append!(idx_mins,m)
        # flatten out data around spikes (and keep indexes in tact)
        if m-ε<1
            dum_df[1:m+ε,2] = I₀*ones(m+ε)
        elseif m+ε>length(dum_df[:,2])
            LD = length(dum_df[:,2])
            dum_df[m-ε:end,2] = I₀*ones(LD-(m-ε)+1)
        else
            dum_df[m-ε:m+ε,2] = I₀*ones(2*ε+1)
        end
    end
    return sort(idx_mins)
end

df2 = deepcopy(df)
idx_min_vec = idx_mins(df2)

## Plot the cropped data
plot_a = plot(df[:,1],df[:,2], title="Light Intensity of a Star", xlabel="Time, Earth Days", ylabel="Intensity", label="Intensity")
scatter!(df[idx_min_vec,1], df[idx_min_vec,2], shape=:star, msize=6, mcolor="red", label="Minima")
savefig(plot_a, "Images/plot_a.png")

## Calculate the orbital period
min_time = df[idx_min_vec,1]
min1 = min_time[1:end-1]
min2 = min_time[2:end]
orb_per = mean(min2 .- min1)


###### Part b: Estimate the ratios of the planetary/solar radii
# Assume Flux from the star is proportional to its area
# Aₚ/Aₛ = (I₀-Iₘ)/I₀

df3 = deepcopy(df)
I₀_vec = reg_intensity(df3)
I₀ = mean(I₀_vec[:,2])
Iₘ = mean(df[idx_min_vec,2])

A_rat = (I₀-Iₘ)/I₀
R_rat = sqrt(A_rat)


###### Part c: Estimate the transit time
## Define functions
function val_drop(mat, δ, skp)
    tran_vec = []
    i = 2
    while i < length(mat[:,2])-1
        if abs(mat[i,2]-mat[i-1,2]) > δ
            append!(tran_vec, i-1)
            i += skp
        else
            i += 1
        end
    end

    return tran_vec
end

function tran_time(mat, δ, skp)
    idx1 = val_drop(mat, δ, skp)
    rev_mat = reverse(mat[1:end-skp,:], dims=1)
    idx2_temp = reverse(val_drop(rev_mat, δ, skp))
    idx2 = length(mat[1:end-skp,1]) .- idx2_temp .+ 1

    L = min(length(idx1), length(idx2))
    idx1 = idx1[1:L]
    idx2 = idx2[1:L]
    time1 = mat[idx1,1]
    time2 = mat[idx2,1]
    time = time2 .- time1

    return mean(time), idx1, idx2
end

## Initialize matrices and execute functions
df4 = deepcopy(df)
mat = hcat(df4[:,1], df4[:,2])
δ = 0.0003*I₀
skp = 10
transit_time, idx1, idx2 = tran_time(mat, δ, skp)

somevec = zeros(Int64, length(idx1))
for i in 1:length(idx1)
    somevec[i] = idx1[i]
end

## Plotting
plot_c = plot(df[:,1],df[:,2], title="Light Intensity of a Star", xlabel="Time, Earth Days", ylabel="Intensity, some units", label="Intensity")
scatter!(df[idx1,1], df[idx1,2], shape=:star, msize=6, mcolor="red", label="Begin Dip")
scatter!(df[idx2,1], df[idx2,2], shape=:star, msize=6, mcolor="green", label="End Dip")
savefig(plot_c, "Images/plot_c.png")
# On my XPS13, for scatter! need to include c=[:red,:red]
