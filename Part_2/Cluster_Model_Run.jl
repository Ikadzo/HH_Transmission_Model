# Author:         Ivy K Kombe
# Institutions:   KEMRI-Wellcome Trust Research Programme, Kilifi, Kenya
#                 London Schoool of Hygiene and Tropical Medicine, London, UK
# Date Published: 23rd September 2019
################################################################################
# This code runs the household transmission model that uses genetic data in the
# form of genetic clusters. The some data used in running the model was modified
# on the R platform and saved as CSV files which a called here. The results of
# the MCMC runs generated here are then saved as CSV files for analysis in R.
################################################################################
# Setting up the environment
############################
# List any variables in the workspace
varinfo()
# Load the functions needed
include("Cluster_Model_Funcs.jl")
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
# For every onset, the start and end day and the person
shed_times_a=Matrix{Int64}(CSV.read("Data_from_R/shed.times.a.csv"))
shed_times_b=Matrix{Int64}(CSV.read("Data_from_R/shed.times.b.csv"))
# The onset time and person index of one case from every HH outbreak
one_case_a1=Matrix{Int64}(CSV.read("Data_from_R/one.case.a.csv"))
one_case_b1=Matrix{Int64}(CSV.read("Data_from_R/one.case.b.csv"))
# Converting into an array of indices
one_case_a=[CartesianIndex(Tuple(one_case_a1[1,:]))]
for i in 2:length(one_case_a1[:,1])
    push!(one_case_a,CartesianIndex(Tuple(one_case_a1[i,:])))
end
one_case_b=[CartesianIndex(Tuple(one_case_b1[1,:]))]
for i in 2:length(one_case_b1[:,1])
    push!(one_case_b,CartesianIndex(Tuple(one_case_b1[i,:])))
end

# ------------ These will change with change in cluster pattern ----------------
# People and dates with sequence information relative to the entire data set
cl_pos_a1=Matrix{Int64}(CSV.read("Data_from_R/cl.pos.a.csv"))
cl_pos_b1=Matrix{Int64}(CSV.read("Data_from_R/cl.pos.b.csv"))
# Converting into an array of indices
cl_pos_a=[CartesianIndex(Tuple(cl_pos_a1[1,:]))]
for i in 2:length(cl_pos_a1[:,1])
    push!(cl_pos_a,CartesianIndex(Tuple(cl_pos_a1[i,:])))
end
cl_pos_b=[CartesianIndex(Tuple(cl_pos_b1[1,:]))]
for i in 2:length(cl_pos_b1[:,1])
    push!(cl_pos_b,CartesianIndex(Tuple(cl_pos_b1[i,:])))
end
# Pair-wise genetic distances between cases with seq data
dist_gen_cases_a=Matrix{Float64}(CSV.read("Data_from_R/dist.gen.cases.a.csv"))
dist_gen_cases_b=Matrix{Float64}(CSV.read("Data_from_R/dist.gen.cases.b.csv"))
# Genetic distance matrices
Nucleotide_D_a=Matrix{Float64}(CSV.read("Data_from_R/Nucleotide.D.a.csv"))
Nucleotide_D_b=Matrix{Float64}(CSV.read("Data_from_R/Nucleotide.D.b.csv"))
# Clustering of sequences
Seq_cluster_id_a=Matrix{Int64}(CSV.read("Data_from_R/Seq.cluster.id.a.csv"))
Seq_cluster_id_b=Matrix{Int64}(CSV.read("Data_from_R/Seq.cluster.id.b.csv"))
# A matrix of shedding episodes identified by cluster type/id for episodes and
# household outbreaks with some genetic info
Shedding_cluster_fixedA=Matrix{Int64}(CSV.read("Data_from_R/Shedding.cluster.fixedA.csv"))
Shedding_cluster_fixedB=Matrix{Int64}(CSV.read("Data_from_R/Shedding.cluster.fixedB.csv"))
# Initial number of HH outbreaks in each cluster
init_rsva_outbreaks=Matrix{Int64}(CSV.read("Data_from_R/init.rsva.outbreaks.csv"))
init_rsvb_outbreaks=Matrix{Int64}(CSV.read("Data_from_R/init.rsvb.outbreaks.csv"))
# Initial random cluster assignment for missing clusters
Shedding_cluster_A=Matrix{Int64}(CSV.read("Data_from_R/Shedding.cluster.a.csv"))
Shedding_cluster_B=Matrix{Int64}(CSV.read("Data_from_R/Shedding.cluster.b.csv"))

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
# - All the people with an onset
shed_pple_a = findall(sum(Onset_a,dims=1).>0);
shed_pple_b = findall(sum(Onset_b,dims=1).>0);
# - Unclustered shedding episodes
x_a = findall((Onset_a .== 1) .& (Shedding_cluster_fixedA .== 0))
x_b = findall((Onset_b .== 1) .& (Shedding_cluster_fixedB .== 0))
# - For every cluster, vectors of unique pair-wise distance
dist_cl_a = Dict([(1,[]), (2,[]), (3,[]), (4,[]), (5,[])])
for cl in 1:maximum(Seq_cluster_id_a)
    dist_cl = Nucleotide_D_a[Seq_cluster_id_a[:] .== cl,Seq_cluster_id_a[:] .== cl];
    if length(dist_cl)>1
        append!(dist_cl_a[cl],unique(dist_cl[findall(UpperTriangular(dist_cl).==0)]));
    else
        append!(dist_cl_a[cl],0);
    end
end
dist_cl_b = Dict([(1,[]), (2,[]), (3,[]), (4,[]), (5,[]), (6,[]),(7,[])])
for cl in 1:maximum(Seq_cluster_id_b)
    dist_cl = Nucleotide_D_b[Seq_cluster_id_b[:] .== cl,Seq_cluster_id_b[:] .== cl];
    if length(dist_cl)>1
        append!(dist_cl_b[cl],unique(dist_cl[findall(UpperTriangular(dist_cl).==0)]));
    else
        append!(dist_cl_b[cl],0);
    end
end
#-------------------------------------------------------------------------------
#                               RUNNING THE MCMC
#-------------------------------------------------------------------------------
# Initial parameter values
theta_names=["PrevHom","PrevHet","CurrHet","age2","age3","age4","etaA","etaB","hhSize","LowSym",
"HighSym","distRate","genRate","epsilonA","epsilonB","epsAge2","epsAge3","delta","beta"]
theta= NamedArray(log.([0.5,0.6,1.0001,1.01,0.3,0.2,0.0188,0.015,0.4,2.48,
                        6.7,100000000,0.0000000000001,0.0033,0.0062,1.0001,2,0.002,0.5]),
                        (theta_names,))
# Given the data and initial parameter values, the log likelihood given by the
# next line is -2129.76
IBM_cluster_LL(theta, Shedding_cluster_A, Shedding_cluster_B)

# --- Chain 1
# A timed test run to check how long it takes to run a few iterations and if the
# algorithm is working
@time(Chain1_runs=MH_MCMC_cluster(Shedding_cluster_A,Shedding_cluster_B,
IBM_cluster_posterior,theta,init_rsva_outbreaks,init_rsvb_outbreaks,25))
# Long runs
Chain1_runs=MH_MCMC_cluster(Shedding_cluster_A,Shedding_cluster_B,
IBM_cluster_posterior,theta,init_rsva_outbreaks,init_rsvb_outbreaks,500000,
NaN,10000,true,1)
# Saving output in CSV format
writedlm( "Res/ParameterTrace1.csv",  Chain1_runs["Parameters"], ',')
writedlm( "Res/TargetTrace1.csv",  Chain1_runs["target.trace"], ',')
writedlm( "Res/ParaAcceptance1.csv",  Chain1_runs["para.acceptance.rate"], ',')
writedlm( "Res/Uninformed_outbreak_a1.csv",  Chain1_runs["uninf.a"], ',')
writedlm( "Res/Uninformed_outbreak_b1.csv",  Chain1_runs["uninf.b"], ',')
writedlm( "Res/AugAccA1.csv",  Chain1_runs["aug.a.acc"], ',')
writedlm( "Res/AugAccB1.csv",  Chain1_runs["aug.b.acc"], ',')

# --- Chain 2
theta_sd=collect(theta)./5
theta2= NamedArray(log.([1.5,1.6,0.5,1.5,1.001,1.2,0.188,0.15,1.4,0.48,
                        0.7,0.0001,1.05,0.033,0.062,0.7,0.2,0.02,0.05]),
                        (theta_names,))
# Long runs
Chain2_runs=MH_MCMC_cluster(Shedding_cluster_A,Shedding_cluster_B,
IBM_cluster_posterior,theta2,init_rsva_outbreaks,init_rsvb_outbreaks,500000,
theta_sd,15000,true,2)
# Saving output in CSV format
writedlm( "Res/ParameterTrace2.csv",  Chain2_runs["Parameters"], ',')
writedlm( "Res/TargetTrace2.csv",  Chain2_runs["target.trace"], ',')
writedlm( "Res/ParaAcceptance2.csv",  Chain2_runs["para.acceptance.rate"], ',')
writedlm( "Res/Uninformed_outbreak_a2.csv",  Chain2_runs["uninf.a"], ',')
writedlm( "Res/Uninformed_outbreak_b2.csv",  Chain2_runs["uninf.b"], ',')
writedlm( "Res/AugAccA2.csv",  Chain2_runs["aug.a.acc"], ',')
writedlm( "Res/AugAccB2.csv",  Chain2_runs["aug.b.acc"], ',')

# --- Chain 3
theta3= theta.+theta2
# Long runs
Chain3_runs=MH_MCMC_cluster(Shedding_cluster_A,Shedding_cluster_B,
IBM_cluster_posterior,theta3,init_rsva_outbreaks,init_rsvb_outbreaks,500000,
theta_sd,10000,true,3)
# Saving output in CSV format
writedlm( "Res/ParameterTrace3.csv",  Chain3_runs["Parameters"], ',')
writedlm( "Res/TargetTrace3.csv",  Chain3_runs["target.trace"], ',')
writedlm( "Res/ParaAcceptance3.csv",  Chain3_runs["para.acceptance.rate"], ',')
writedlm( "Res/Uninformed_outbreak_a3.csv",  Chain3_runs["uninf.a"], ',')
writedlm( "Res/Uninformed_outbreak_b3.csv",  Chain3_runs["uninf.b"], ',')
writedlm( "Res/AugAccA3.csv",  Chain3_runs["aug.a.acc"], ',')
writedlm( "Res/AugAccB3.csv",  Chain3_runs["aug.b.acc"], ',')
