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
function dum(df)
    idx_mins = []
    ε = 5
    dum_df = copy(df)
    for i in 1:5 # number of dips (obvious from graph)
        m = argmin(dum_df[:,2])
        append!(idx_mins,m)
        dum_df = dum_df[setdiff(1:end, m-ε:m+ε), :] # remove the rows of data around the minimum
    end
    return dum_df, sort(idx_mins)
end

I₀_vec, idx_min_vec = dum(df)

## Plot the cropped data
plot(I₀_vec[:,1],I₀_vec[:,2], title="Light Intensity of a Star", xlabel="Time, Earth Days", ylabel="Intensity, some units",legend=false)

## Calculate the orbital period
min_time = df[idx_min_vec,1]
min1 = min_time[1:end-1]
min2 = min_time[2:end]
orb_per = mean(min2 .- min1)


###### Part b: Estimate the ratios of the planetary/solar radii
# Assume Flux from the star is proportional to its area
# Aₚ/Aₛ = (I₀-Iₘ)/I₀

I₀ = mean(I₀_vec[:,2])
Iₘ = mean(df[idx_min_vec,2])

rat = (I₀-Iₘ)/I₀


###### Part c: Estimate the transit time
function val_drop(mat, δ, skp)
    tran_vec = []
    i = 2
    while i < length(mat[:,2])-1
        println(i)
        if abs(mat[i,2]-mat[i-1,2]) > δ |
            println("skip at i = $i")
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

    return mean(time)
end

mat = hcat(df[:,1], df[:,2])
δ = 0.0003*I₀
skp = 10
transit_time = tran_time(mat, δ, skp)

plot(df[:,1],df[:,2], title="Light Intensity of a Star", xlabel="Time, Earth Days", ylabel="Intensity, some units",legend=false)
scatter!(df[idx1,1], df[idx1,2], shape=:star, msize=8, mcolor="red")
scatter!(df[idx2,1], df[idx2,2], shape=:star, msize=8, mcolor="green")
