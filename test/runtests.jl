Pkg.activate(".")
import SINDYPI_Julia

using SINDYPI_Julia

using Test

#Test 1: Glycolysis
using MAT
vars = matread("C:/Users/Raj/Desktop/SINDY PI Brunton/SINDy-PI-master/SINDy-PI-master/Comparison/DataLength/YeastGlycolysis/SINDy_PI/Datas/TrainingData.mat")
X = vars["xt"]
DX = vars["dxt"]

Polyorder =  Int64[6 6 2 2 2 6 2]
Trigorder = Int64[0 0 0 0 0 0 0]

SINDYPI(X, DX, Polyorder, Trigorder)

#Test 2: Pendulum on Cart
using DifferentialEquations, Plots
# Setup
initial = [0.3, 0, 1.0, 0]
tspan = (0,16.0)
dt = 0.001;

 m = 1
 M = 1
 L = 1
 g = 9.81


# Define the function
function Cart_Pendulum!(du, u,p,t)
    du[1] = u[3]
    du[2] = u[4]
    du[3] = -((M +m)*g*sin(u[1]) +m*(L^2)*sin(u[1])*cos(u[1])*u[3]*u[3] )/(L*L*(M + m -m*cos(u[1])*cos(u[1])))
    du[4] = (m*L*L*sin(u[1])*u[3]*u[3] + m*g*sin(u[1])*cos(u[1]))/(L*(M + m -m*cos(u[1])*cos(u[1])))
    return [du[1] du[2] du[3] du[4]]
end

#Pass to solvers
prob = ODEProblem(Cart_Pendulum!, initial, tspan)

sol = solve(prob, Tsit5(), saveat = dt)

plot(sol,linewidth=2,
     title ="Pendulum on Cart Problem",
     xaxis = "Time",
     yaxis = "Height",
     label = ["Theta","S", "dTheta", "dS"])

X = Array(sol)
DX = similar(X)
du = zeros(4,1)
for (i, xi) in enumerate(eachcol(X))
    DX[:,i] = Cart_Pendulum!(du, xi, [], 0.0)
end
X = X'
DX = DX'

Polyorder =  Int64[1 1 0 0]
Trigorder = Int64[0 0 3 4]

SINDYPI(X, DX, Polyorder, Trigorder)

#TEST 3: Michaelis - Menten
jx = 0.6
Vm = 1.5
Km = 0.3

f(u,p,t) = jx - (Vm*u/(Km + u))

u0 = 1/2
tspan = (0.0,3.0)
dt = 0.01
prob = ODEProblem(f,u0,tspan)

sol = DifferentialEquations.solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8, saveat = dt)

X = sol.u
DX = jx .- (Vm*X./(Km .+ X))

Polyorder =  Int64[1]
Trigorder = Int64[0]

SINDYPI(X, DX, Polyorder, Trigorder)




@testset "SINDYPI_Julia.jl" begin
    # Write your own tests here.
end
