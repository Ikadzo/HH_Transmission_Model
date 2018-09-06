Project: Model-based estimates of transmission of respiratory syncytial virus within households.

This work describes the analysis of household cohort data collected over a period of 6 months covering an epidemic of respiratory syncytial virus (RSV). The data was collected from a rural coastal community in Kenya in 2009.The aim of the analysis is to gain a better understanding of RSV transmission for the purpose of intervention planning. We show that household size, symptom status and viral load are important factors of within household transmission. We find evidence of the possible existence of an RSV group specific transmission niche that could form part of the explanation for RSV A and B temporal and geographical co-existence.  

General instructions
These instructions will allow you to be able to replicate the results shown in the main research article. 

Prerequisites
The analysis was done in RStudio Version 1.1.383 - https://www.rstudio.com/. The following packages are required: 


Dataset 
The dataset is available from ... 

File list
The following are the files that accompany the data:
	1. Code scripts used to run the analysis and generate results
		a) Data_Prep.R - 
		b) Model_Run.R - 
		c) Model_Diag.R - 
		d) Model_Funcs.R
		e) Sampler_fitR.R - 
	2. Variable codebook - describing variables contained in the dataset including data types and value labels

Running the code
Step 1:
Download the files into your designated folder

Step 2:
Run the main code 'Model_Run.R'. This function has two parts: PART 1 runs the model while distinguishing between RSV A and B, PART 2 does not. If you wish to only run one part, place break points or comment out the relevant sections of the code.Each part runs 3 MCMC chains sequentially. PART 1 is set to run 250,000 iterations while PART 2 is set to run 100,000 but these can be changed in the code. The diagnostic process is interactive and will require the user to input a burn-in point and thinning interval for each chain.


