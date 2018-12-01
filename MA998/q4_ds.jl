#################################################
################# Question 4 ####################
#################################################

## Packages
using DynamicalSystems

## Use Julia Plots, with PyPlot backend
using Plots
pyplot()

## define maps as a dynamical system
function g(dx, x, p, n)
    a, b, c, d = p
    dx[1] = x[1]^2 - x[2]^2 + a*x[1] + b*x[2]
    dx[2] = 2*x[1]*x[2] + c*x[1] + d*x[2]
    return dx
end

function g_jac(J, x, p, n)
    J[1,1] = 2*x[1]+p[1]
    J[1,2] = -2*x[2]+p[2]
    J[2,1] = 2*x[2]+p[3]
    J[2,2] = 2*x[1]+p[4]
    return
end

st = sqrt(2)/10
p_1 = [0.9, -0.6013, 2, 0.5]
p_2 = [0.3, 0.6, 2, 0.27]
ds_1 = DiscreteDynamicalSystem(g, [st, st], p_1, g_jac)
ds_2 = DiscreteDynamicalSystem(g, [st, st], p_2, g_jac)

## get trajectory of the dynamical system
traj_1 = trajectory(ds_1, 10000)
x_1 = traj_1[:,1]
y_1 = traj_1[:,2]

traj_2 = trajectory(ds_2, 10000)
x_2 = traj_2[:,1]
y_2 = traj_2[:,2]

## plotting
chaos = scatter(x_1,y_1, msize=1, title="First Set of Parameters", xlabel="x_n", ylabel="y_n", legend=false)
savefig(chaos, "Images/chaos.png")
ring = scatter(x_2,y_2, msize=1, title="Second Set of Parameters", xlabel="x_n", ylabel="y_n", legend=false)
savefig(ring, "Images/ring.png")
