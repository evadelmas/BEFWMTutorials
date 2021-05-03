#=
This is the accompanying script for the Rewiring in the Bioenergetic Food-Web Model (BEFWM)
notebook
=#

import Pkg
Pkg.activate(".")
using BioEnergeticFoodWebs, EcologicalNetworks
include("src/utils.jl") 
import Random.seed!
import Statistics.mean
seed!(55586)
Pkg.status()

# Step 1

S = 30 #number of species
C = 0.15 #expected connectance
A = EcologicalNetworks.nichemodel(S, C).edges #generates the food web
A = Int.(Array(A)) #transform the interaction matrix into a format compatible with BEFWM
pltA = webplot(A, consasrow = true) #true means that we want to have consumers plotted as rows(i) and resources as columns (j)

S = 30 #number of species
C = 0.15 #expected connectance
A = EcologicalNetworks.nichemodel(S, C).edges #generates the food web
A = Int.(Array(A)) #transform the interaction matrix into a format compatible with BEFWM
pltA = webplot(A, consasrow = true) #true means that we want to have consumers plotted as rows(i) and resources as columns (j)

# Generate the sets of parameters

#no rewiring is the default, technically you don't need to specify rewire_method in that case
p_none = model_parameters(A, bodymass = M, h = 2.0, rewire_method = :none) 

# Set up simulations
b0 = rand(S) #random starting biomasses
tstop = 10000 #simulation time
ϵ = 1e-10 #extinction threshold
#when a species biomass reaches the extinction threshold, it's considered as
#extinct. It's important to know the extinction threshold when using rewiring because 
#rewiring will be triggered at each extinction event during the simulations. 
# (We'll see later how this can be different for the ADBM)

# Step 2

## Diet Overlap

p_do = model_parameters(A, bodymass = M, h = 2.0, rewire_method = :DO) 
s_do = simulate(p_do, b0, stop = tstop, extinction_threshold = ϵ)

plt_dyn_do = plot(s_do[:t], s_do[:B], leg = false, c = :black, ylims = (-0.01,3), xlabel = "time", ylabel = "biomass")
plt_mat_do = webplot(p_do[:A], consasrow = true)
plot(plt_dyn_do, plt_mat_do, size = (750, 400))

## Diet SImilarity 

#TODO

## Allometric Diet Breadth

#TODO

# Step 3

# This function (defined in utils.jl) does exactly what is axplained above, 
# it removes the species that are extinct from the interaction matrix
# and also checks for disconnected consumers. If this function detects any disconnected consumer
# it will return a vector of species identity (the disconnected consumers identified) instead of 
# an updated interaction matrix
new_A = updateA(s_do)
webplot(new_A)

original_height = maximum(trophic_rank(A))
updated_height = maximum(trophic_rank(new_A))
original_avgindegree = mean(sum(A, dims = 1))
updated_avgindegree = mean(sum(new_A, dims = 1))
original_connectance = sum(A) / (S ^ 2)
updated_connectance = sum(new_A) / (size(new_A, 1) ^ 2)
structural_change = (height = round(updated_height - original_height, digits = 2)
                   , avg_indegree = round(updated_avgindegree - original_avgindegree, digits = 2)
                   , connectance = round(updated_connectance - original_connectance, digits = 2)) 