# Author:         Ivy K Kombe
# Institutions:   KEMRI-Wellcome Trust Research Programme, Kilifi, Kenya
#                 London Schoool of Hygiene and Tropical Medicine, London, UK
# Date Published: 24th September 2019
################################################################################
# This code contains functions needed for running the model that uses data
# identified at the pathogen level for the purpose of making comparisons
# to the model in 'Cluster_Model_Run.jl'. There are 5 functions in
################################################################################
# Setting up the environment
############################
# List variables in the workspace
varinfo()
# Load the functions needed
include("Pathogen_Model_Funcs.jl")
# Load the libraries needed
using CSV; using LinearAlgebra; using DataFrames;using NamedArrays;using DelimitedFiles
using StatsBase;using Distributions;

###################################
# Reading in the data modified in R
###################################
# 0/1 matrix of shedding episodes for RSV A and B
Shedding_rsv=Matrix{Int64}(CSV.read("Data_from_R/Shedding.rsv.csv"));
# Categorical variable for infection history (3 groups)
B=Matrix{Int64}(CSV.read("Data_from_R/B.rsv.csv"))
B[B .== 0].=1;
# Combined categories of viral load and symptoms
Load_n_sym_rsv=Matrix{Int64}(CSV.read("Data_from_R/Load.n.sym.rsv.csv"))
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
Time=length(Shedding_rsv[:,1])
Population=length(Shedding_rsv[1,:]);
# - Matrix of onset days given shedding matrix
# RSV
Onset_rsv=zeros(Time-1,Population);x=findall(isequal(1),diff(Shedding_rsv,dims=1));
Onset_rsv[x] .= 1;Onset_rsv=vcat(0,Onset_rsv)
# - Matrix of at risk status, i.e. present and not shedding
atRisk_rsv = copy(status); atRisk_rsv[Shedding_rsv .== 1]= repeat([0],sum(Shedding_rsv .== 1));
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
Shedding_cluster=Shedding_rsv;

# Initial parameter values
theta_names=["PrevHom","age2","age3","age4","eta","hhSize","LowSym",
"HighSym","distRate","epsilon","epsAge2","epsAge3","delta","beta"]
theta= NamedArray(log.([0.5,1.01,0.3,0.2,0.0188,0.4,2.48,
                        6.7,1000000,0.0033,1.0001,2,0.002,0.5]),
                        (theta_names,))
# --- Chain 1
# Long runs
Chain1_runs=MH_MCMC(IBM_posterior,theta,150000,NaN,10000,true,1)
# Saving output in CSV format
writedlm( "Res/ParameterTraceRSV1.csv",  Chain1_runs["Parameters"], ',')
writedlm( "Res/TargetTraceRSV1.csv",  Chain1_runs["target.trace"], ',')
writedlm( "Res/ParaAcceptanceRSV1.csv",  Chain1_runs["para.acceptance.rate"], ',')


# Chain 2
#----------
theta_sd=collect(theta)./5
theta2= NamedArray(log.([1.5,2,1.3,1.2,0.00188,1.4,0.48,
                        0.7,1.000001,0.033,1.6,0.2,0.02,0.05]),
                        (theta_names,))
Chain2_runs=MH_MCMC(IBM_posterior,theta2,150000,theta_sd,10000,true,2)

writedlm( "Res/ParameterTraceRSV2.csv",  Chain2_runs["Parameters"], ',')
writedlm( "Res/TargetTraceRSV2.csv",  Chain2_runs["target.trace"], ',')
writedlm( "Res/ParaAcceptanceRSV2.csv",  Chain2_runs["para.acceptance.rate"], ',')

# Chain 3
#----------
theta3= theta.+theta2
Chain3_runs=MH_MCMC(IBM_posterior,theta3,150000,theta_sd,10000,true,3)

writedlm( "Res/ParameterTraceRSV3.csv",  Chain3_runs["Parameters"], ',')
writedlm( "Res/TargetTraceRSV3.csv",  Chain3_runs["target.trace"], ',')
writedlm( "Res/ParaAcceptanceRSV3.csv",  Chain3_runs["para.acceptance.rate"], ',')
