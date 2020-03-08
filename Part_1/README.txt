Author:         Ivy K Kombe
Institutions:   KEMRI-Wellcome Trust Research Programme, Kilifi, Kenya
                London School of Hygiene and Tropical Medicine, London, UK
Date Published: 13th September 2018

Project:        Model-based estimates of transmission of respiratory syncytial virus within households.

This work describes the analysis of household cohort data collected over a period of 6 months covering an epidemic of respiratory syncytial virus (RSV). The data was collected from a rural coastal community in Kenya in 2009.The aim of the analysis is to gain a better understanding of RSV transmission for the purpose of intervention planning. We show that household size, symptom status and viral load are important factors of within household transmission. We find evidence of the possible existence of an RSV group specific transmission niche that could form part of the explanation for RSV A and B temporal and geographical co-existence.  

General instructions
These instructions will allow you to be able to replicate the results shown in the main research article. 

Prerequisites
The analysis was done in RStudio Version 1.1.383 - https://www.rstudio.com/. The following packages are required: 


Dataset 
The data is available under the Creative Commons Attribution 4.0 International (CC BY 4.0) - 
https://creativecommons.org/licenses/by/4.0/legalcode. For access to data and more detailed information beyond the metadata provided, there is a process of managed access requiring submission of a request form for consideration by the Data Governance Committee at KEMRI-Wellcome Trust Research Programme (http://kemri-wellcome.org/about-us/#ChildVerticalTab_15). There are two data files provided for the analysis. One is the household cohort data. The other is on the results of RSV antigen testing conducted at the Kilifi County Hospital during the same time period as the household study. These data are for children <5 years of age. 

File list
The following are the files that accompany the data:
	1. Code scripts used to run the analysis and generate results
		a) Data_Prep.R - Modifies the data for fitting
		b) Model_Run.R - Runs the MCMC fitting algorithm
		c) Model_Diag.R - Runs diagnostics on the results of the MCMC chains
		d) Model_Funcs.R - Contains the prior, likelihood and posterior functions
		e) Sampler_fitR.R - The package 'fitR' that contains the MCMC functions needed will not load on the latest version of R. As a work-around, I copied the function codes from github.
		f) Model_Validate - Runs the model validation process
		g) Model_SA - Runs the sensitivity analysis described in the main article

	2. Variable codebook - describing variables contained in the dataset including data types and value labels

Running the main code
Step 1:
Download the files into your designated folder

Step 2:
Run the main code 'Model_Run.R'. This function has two parts: PART 1 runs the model while distinguishing between RSV A and B, PART 2 does not. If you wish to only run one part, place break points or comment out the relevant sections of the code.Each part runs 3 MCMC chains sequentially. PART 1 is set to run 250,000 iterations while PART 2 is set to run 100,000 but these can be changed in the code. The diagnostic process is interactive and will require the user to input a burn-in point and thinning interval for each chain.


