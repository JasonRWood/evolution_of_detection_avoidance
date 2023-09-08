import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp
import math


sys.path.append("../src/Model_B_gillespie")

import runner_gillespie_B

solB = runner_gillespie_B.PySolver()

N = 100000

starting_infecteds = 10

resolution = 10

# zetas = [0, 3*5, 3*10, 3*15, 3*20]
param_res = 50
zeta_max = 50
zetas = [zeta_max*i/param_res for i in range(param_res + 1)]

beta_max = 0.0001
c1 = 1.0
c2 = -3

S_increment = 0


repeats = min(mp.cpu_count(), 8)

t_max = 100

eta = 0.9

delta = 0.1

alpha = 0.5
gamma = 0.5

seed = 5000
np.random.seed(seed)

mut_chance = 0.2

simulated_rhos = []

seeds = []
for i in range(len(zetas)):
    temp_seeds = []
    for j in range(repeats):
        temp_seeds.append(seed + j + i*repeats)
    seeds.append(temp_seeds)
    
# print(seeds[-1][-1])
# flop
average_rhos = {}
num_infecteds_average = {}
t_step_counts = {}
plotting_rhos = []
plotting_zetas = []
data = []
for i, zeta in enumerate(zetas):
    average_rho = []
    num_infecteds_average[zeta] = {}
    
    t_step_counts[zeta] = {}
    processes = []
    temp_rhos = []
    for j in range(repeats):
        p = mp.Process(target = solB.gillespie_simulation, args = (seeds[i][j], N, starting_infecteds, beta_max, c1, c2, zeta, eta, delta, alpha, gamma, mut_chance, t_max, resolution, S_increment))
        p.start()
        processes.append(p)
    for proc in processes:
        proc.join()
        solB.gillespie_simulation(seed, N, starting_infecteds, beta_max, c1, c2, zeta, eta, delta, alpha, gamma, mut_chance, t_max, resolution, S_increment)
    for j in range(repeats):
        df = pd.read_csv(f"../data/gillespie_simulations_B/data_set{seeds[i][j]}.csv")
        df = df[df["t_step"] <= t_max]
        
#         seed += 1
        
        evolved_rhos = [val for val in df["Average_rho"].values]
        t_steps = [val for val in df["t_step"].values]
        infecteds = [val for val in df["num_infected"].values]
        
        max_infecteds = max(infecteds)
        if max_infecteds < 1000:
            print(seeds[i][j], zeta, max_infecteds)
        else:
            inf_ind = infecteds.index(max_infecteds)
            rho_max = evolved_rhos[inf_ind]
            plotting_rhos.append(rho_max)
            temp_rhos.append(rho_max)
            plotting_zetas.append(zeta)
    data.append(temp_rhos)
    
fig, ax = plt.subplots(figsize=(10, 5))
 
# Creating plot
bp = ax.boxplot(data)

# plt.scatter(plotting_zetas, plotting_rhos)
desired_ticks = [int((i/5)*param_res) for i in range(6)]
labels = []
ticks = []
for i in range(param_res+1):
    ticks.append(i + 1)
    if i in desired_ticks:
        labels.append(zetas[i])
    else:
        labels.append("")
plt.xticks(ticks, labels)
plt.ylim(bottom = 0)
plt.xlabel(r"Testing rate, $\zeta$")
plt.ylabel(r"Average detection avoidance, $\rho$")
plt.savefig("../outputs/Model_B_parameter_sweep_bp.png", bbox_inches = "tight")
plt.savefig("../outputs/Model_B_parameter_sweep_bp.pdf", bbox_inches = "tight")
plt.close()