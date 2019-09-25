# Author:         Ivy K Kombe
# Institutions:   KEMRI-Wellcome Trust Research Programme, Kilifi, Kenya
#                 London Schoool of Hygiene and Tropical Medicine, London, UK
# Date Published: 24th September 2019
################################################################################
# This code contains functions needed for running the model that uses data
# identified at the pathogen group level for the purpose of making comparisons
# to the model in 'Cluster_Model_Run.jl'. Unlike the model in the publication
# 'Model based estimates of transmission of respiratory syncytial virus within
# households' this version modifies susceptibility based on current shedding
# status and uses an emperical based background community function.
# Some of the data used in running the model was modified
# on the R platform and saved as CSV files which are called here. The results of
# the MCMC runs generated here are then saved as CSV files for analysis in R.
################################################################################
# Setting up the environment
############################
# List variables in the workspace
varinfo()
# Load the functions needed
include("Group_Model_Funcs.jl")
# Load the libraries needed
using CSV; using LinearAlgebra; using DataFrames;using NamedArrays;using DelimitedFiles
using StatsBase;using Distributions;

###################################
# Reading in the data modified in R
###################################
# 0/1 matrix of shedding episodes for RSV A and B
Shedding_a=Matrix{Int64}(CSV.read("Data_from_R/Shedding.a.csv"));
Shedding_b=Matrix{Int64}(CSV.read("Data_from_R/Shedding.b.csv"));
# Categorical variable for infection history (3 groups)
B=Matrix{Int64}(CSV.read("Data_from_R/B.csv"))
# Combined categories of viral load and symptoms
Load_n_sym_a=Matrix{Int64}(CSV.read("Data_from_R/Load.n.sym.a.csv"))
Load_n_sym_b=Matrix{Int64}(CSV.read("Data_from_R/Load.n.sym.b.csv"))
# 0/1 matrix of presence or absence in the household
status=Matrix{Int64}(CSV.read("Data_from_R/status.csv"))
# Demographic information
Demo_data=CSV.read("Data_from_R/Demo.data.csv")
# Spatial distance matrix
Distance_matrix=Matrix{Float64}(CSV.read("Data_from_R/Distance.matrix.csv"))

#########################################################
# Genarating other variables needed for running the model
#########################################################
# - Latency distribution
Lat_dist = DataFrame(days=[0,1,2,3,4,5],prob=[0, 0, 4, 4, 3, 1]/12)
# - Size of the data
Time=length(Shedding_a[:,1])
Population=length(Shedding_a[1,:]);
# - Matrix of onset days given shedding matrix
# RSV A
Onset_a=zeros(Time-1,Population);x=findall(isequal(1),diff(Shedding_a,dims=1));
Onset_a[x] .= 1;Onset_a=vcat(0,Onset_a)
# RSV B
Onset_b=zeros(Time-1,Population);x=findall(isequal(1),diff(Shedding_b,dims=1));
Onset_b[x].=1;Onset_b=vcat(0,Onset_b)
# - Matrix of at risk status, i.e. present and not shedding
atRisk_a = copy(status); atRisk_a[Shedding_a .== 1]= repeat([0],sum(Shedding_a .== 1));
atRisk_b = copy(status); atRisk_b[Shedding_b .== 1]= repeat([0],sum(Shedding_b .== 1));
# - An NxN matrix the identifies housemates.Used in FOI function to add up HH
#   shedders
houses=unique(Demo_data.hhid);
Cc=zeros(Population,Population)
for n in 1:length(houses)
    # Their housemates
    pp=findall(isequal(houses[n]),Demo_data.hhid)
    Cc[pp,pp]=repeat([1],length(pp)^2)
end
# replacing the diagonal with 0's
Cc[diagind(Cc)]=repeat([0],Population)
# - An NxN matrix the identifies non-housemates.Used in FOI function to add up
#   neighbour shedders
Dd = 0 .* Cc;Dd[Cc .== 0] = repeat([1],sum(Cc .== 0))
Dd[diagind(Dd)] =repeat([0],Population)
# - HH size covariate
H = repeat([2],Population);H[Demo_data.hhsize .< 8] .= 1;
# -  Age covariate for susceptibility
A = repeat([4],Population);
A[Demo_data.ageyrs .< 1].=1;
A[(Demo_data.ageyrs.>=1) .& (Demo_data.ageyrs.<5)].= 2;
A[(Demo_data.ageyrs.>=5) .& (Demo_data.ageyrs.<15)].= 3;
# - Age Covariate for community risk/exposure
E=repeat([3],Population);
E[Demo_data.ageyrs .< 1] .= 1;
E[(Demo_data.ageyrs.>=1) .& (Demo_data.ageyrs.<5)] .= 2;
#-------------------------------------------------------------------------------
#                               RUNNING THE MCMC
#-------------------------------------------------------------------------------
Shedding_cluster_A=Shedding_a;
Shedding_cluster_B=Shedding_b;
# Initial parameter values
theta_names=["PrevHom","PrevHet","CurrHet","age2","age3","age4","etaA","etaB","hhSize","LowSym",
"HighSym","distRate","epsilonA","epsilonB","epsAge2","epsAge3","delta","beta"]
theta= NamedArray(log.([0.5,0.6,1.0000001,1.01,0.3,0.2,0.0188,0.015,0.4,2.48,
                        6.7,1000000,0.0033,0.0062,1.0001,2,0.002,0.5]),
                        (theta_names,))

# --- Chain 1
# Long runs
Chain1_runs=MH_MCMC(IBM_posterior,theta,250000,NaN,10000,true,1)
# Saving output in CSV format
writedlm( "Res/ParameterTraceNull1.csv",  Chain1_runs["Parameters"], ',')
writedlm( "Res/TargetTraceNull1.csv",  Chain1_runs["target.trace"], ',')
writedlm( "Res/ParaAcceptanceNull1.csv",  Chain1_runs["para.acceptance.rate"], ',')


# --- Chain 2
theta_sd=collect(theta)./5
theta2= NamedArray(log.([1.5,1.6,0.5,1.5,1.001,1.2,0.188,0.15,1.4,0.48,
                        0.7,0.0001,0.033,0.062,0.7,0.2,0.02,0.05]),
                        (theta_names,))
# Long runs
Chain2_runs=MH_MCMC(IBM_posterior,theta2,250000,theta_sd,10000,true,2)
# Saving output in CSV format
writedlm( "Res/ParameterTraceNull2.csv",  Chain2_runs["Parameters"], ',')
writedlm( "Res/TargetTraceNull2.csv",  Chain2_runs["target.trace"], ',')
writedlm( "Res/ParaAcceptanceNull2.csv",  Chain2_runs["para.acceptance.rate"], ',')

# --- Chain 3
theta3= theta.+theta2
# Long runs
Chain3_runs=MH_MCMC(IBM_posterior,theta3,250000,theta_sd,10000,true,3)
# Saving output in CSV format
writedlm( "Res/ParameterTraceNull3.csv",  Chain3_runs["Parameters"], ',')
writedlm( "Res/TargetTraceNull3.csv",  Chain3_runs["target.trace"], ',')
writedlm( "Res/ParaAcceptanceNull3.csv",  Chain3_runs["para.acceptance.rate"], ',')
