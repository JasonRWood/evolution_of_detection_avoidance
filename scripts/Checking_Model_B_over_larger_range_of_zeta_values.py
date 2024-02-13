"""
This script uses the python multiprocessing library to perform multiple stochastic simulations
of the large non-replenishing population. The average level of detection avoidance at the peak number
of infected individuals is then observed. After which box-plots are created to observe the variance
of these simulations
"""

# Importing necessary libraries
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp
import math

#Reading in the C++ wrapper file
sys.path.append("../src/Model_B_gillespie")

import runner_gillespie_B

solB = runner_gillespie_B.PySolver()

#Setting our population and setting the number of infecteds
N = 100000

starting_infecteds = 10

#Phenotype space resolution
resolution = 10

# How many zeta values we check, and what the maximum zeta value is 
param_res = 50
zeta_max = 50
zetas = [zeta_max*i/param_res for i in range(param_res + 1)]

#Parameters of the trade-off
beta_max = 0.0001
c1 = 1.0
c2 = -3

# S_increment is used to determine whether or not the population replenishes itself
S_increment = 0

#The number of repeats
repeats = 8

#The number of available cores
num_cpu = mp.cpu_count()

# Other simulation parameters
t_max = 100

eta = 0.9

delta = 0.1

alpha = 0.5
gamma = 0.5

#Seed base number
seed = 5000
np.random.seed(seed)

mut_chance = 0.2

simulated_rhos = []

#Here we pre-generate the seeds we will need, due to a quirk of the multiprocessing library
seeds = []
for i in range(len(zetas)):
    temp_seeds = []
    for j in range(repeats):
        temp_seeds.append(seed + j + i*repeats)
    seeds.append(temp_seeds)
    


#Dictonaries and arrays used to store outputs from simulations
average_rhos = {}
num_infecteds_average = {}
t_step_counts = {}
plotting_rhos = []
plotting_zetas = []
data = []

#Performing the simulations for each zeta value
for i, zeta in enumerate(zetas):
    average_rho = []
    num_infecteds_average[zeta] = {}
    
    t_step_counts[zeta] = {}
    processes = []
    temp_rhos = []
    cpu_count = 0
    #Looping through our repeats and initialising a simulation for each value
    for j in range(repeats):
        if cpu_count <= num_cpu:
            p = mp.Process(target = solB.gillespie_simulation, args = (seeds[i][j], N, starting_infecteds, beta_max, c1, c2, zeta, eta, delta, alpha, gamma, mut_chance, t_max, resolution, S_increment))
            p.start()
            processes.append(p)
            cpu_count += 1
        else:
            for proc in processes:
                proc.join()
            p = mp.Process(target = solB.gillespie_simulation, args = (seeds[i][j], N, starting_infecteds, beta_max, c1, c2, zeta, eta, delta, alpha, gamma, mut_chance, t_max, resolution, S_increment))
            p.start()
            processes = [p]
            cpu_count = 1
    
    #Collecting the processes to ensure all are finished before we read in the csvs
    for proc in processes:
        proc.join()
        
    #Read in and analyse the data from the simulation
    for j in range(repeats):
        
        #Reading in data and filtering
        df = pd.read_csv(f"../data/gillespie_simulations_B/data_set{seeds[i][j]}.csv")
        df = df[df["t_step"] <= t_max]
        
        #Making data accesible in lists for personal preference
        evolved_rhos = [val for val in df["Average_rho"].values]
        t_steps = [val for val in df["t_step"].values]
        infecteds = [val for val in df["num_infected"].values]
        
        #Calculating the number of maximum infecteds
        max_infecteds = max(infecteds)
        
        #Discarding simulations where our threshold was not met
        if max_infecteds < 1000:
            print(seeds[i][j], zeta, max_infecteds)
        #Otherwise determining the level of detection avoidance at the infection peak
        else:
            inf_ind = infecteds.index(max_infecteds)
            rho_max = evolved_rhos[inf_ind]
            plotting_rhos.append(rho_max)
            temp_rhos.append(rho_max)
            plotting_zetas.append(zeta)
            
    #Adding our data
    data.append(temp_rhos)
    
#Creating Figure
fig, ax = plt.subplots(figsize=(10, 5))
 
# Creating plot
bp = ax.boxplot(data)

# Making appropriate ticks
desired_ticks = [int((i/5)*param_res) for i in range(6)]
labels = []
ticks = []
for i in range(param_res+1):
    ticks.append(i + 1)
    if i in desired_ticks:
        labels.append(zetas[i])
    else:
        labels.append("")
        
#Updating ticks and saving figure
plt.xticks(ticks, labels)
plt.ylim(bottom = 0)
plt.xlabel(r"Testing rate, $\zeta$")
plt.ylabel(r"Average detection avoidance, $\rho$")
plt.savefig("../outputs/Model_B_parameter_sweep_bp.png", bbox_inches = "tight")
plt.savefig("../outputs/Model_B_parameter_sweep_bp.pdf", bbox_inches = "tight")
plt.close()