import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
import math

sys.path.append("../src/Model_A_gillespie")

import runner_gillespie_A

sys.path.append("../src/Model_B_gillespie")

import runner_gillespie_B

sys.path.append("../src/Model_C_gillespie")

import runner_gillespie_C

solA = runner_gillespie_A.PySolver()

solB = runner_gillespie_B.PySolver()

solC = runner_gillespie_C.PySolver()

N = 100000

starting_infecteds = 10

resolution = 10

zetas = [0, 3*1, 3*2, 3*3, 3*4]
sigmas = [0.0, 0.25, 0.5, 0.75, 1.0]

labels = [["(A)", "(D)"], ["(B)", "(E)"], ["(C)", "(F)"]]

beta_max = 0.0001
c1 = 1.0
c2 = -3

S_increment = 0

repeats = 20

t_max = 100

eta = 0.9

delta = 0.1

alpha = 0.5
gamma = 0.5

seed = 5000
np.random.seed(seed)

mut_chance = 0.2

colours = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:gray"]

steps = 100000
    
fig = plt.figure(figsize = (15, 10))
gs = fig.add_gridspec(2,3)
ax = [[fig.add_subplot(gs[i,j]) for i in range(2)] for j in range(3)]

recording_size = 1000
recording_mesh = [i/recording_size for i in range(recording_size + 1)]
# print(recording_mesh)
# flop
print("Doing A")
average_rhos = {}
num_infecteds_average = {}
t_step_counts = {}
for i, sigma in enumerate(sigmas):
    average_rhos[sigma] = {}
    num_infecteds_average[sigma] = {}
    t_step_counts[sigma] = {}
    num_sims = 0
    for j in range(repeats):
        solA.gillespie_simulation(seed, N, starting_infecteds, beta_max, c1, c2, sigma, eta, delta, alpha, gamma, mut_chance, t_max, resolution, S_increment)

        df = pd.read_csv(f"../data/gillespie_simulations_A/data_set{seed}.csv")
        
        seed += 1
        
        df = df[df["t_step"] <= t_max]
        
        evolved_rhos = [val for val in df["Average_rho"].values]
        t_steps = [val for val in df["t_step"].values]
        infecteds = [val for val in df["num_infected"].values]
        
        max_infecteds = max(infecteds)
        if max_infecteds > 1000:
            num_sims += 1
            inf_ind = infecteds.index(max_infecteds)
            t_peak = t_steps[inf_ind]
            recast_ts = [t/t_peak for t in t_steps[:inf_ind+1]]
    #         print(recast_ts[-1], t_peak)
            recast_infecteds = []
            recast_rhos = []
            rec_ind = 0
            for k, val in enumerate(recast_ts):
                if val >= recording_mesh[rec_ind]:
                    try:
                        average_rhos[sigma][recording_mesh[rec_ind]] += evolved_rhos[k]
                    except KeyError:
                        average_rhos[sigma][recording_mesh[rec_ind]] = evolved_rhos[k]

                    try:
                        num_infecteds_average[sigma][recording_mesh[rec_ind]] += infecteds[k]
                    except KeyError:
                        num_infecteds_average[sigma][recording_mesh[rec_ind]] = infecteds[k]

                    try:
                        t_step_counts[sigma][recording_mesh[rec_ind]] += 1
                    except KeyError:
                        t_step_counts[sigma][recording_mesh[rec_ind]] = 1


                    rec_ind += 1
        
#         print(t_peak)
#         flop
#         for k, val in enumerate(evolved_rhos):
#             if val > 1:
#                 print(seed, t_steps[k], val, infecteds[k])
#         for k, val in enumerate(t_steps):
#             try:
#                 average_rhos[sigma][val] += evolved_rhos[k]
#             except:
#                 average_rhos[sigma][val] = evolved_rhos[k]
            
#             try:
#                 num_infecteds_average[sigma][val] += infecteds[k]
#             except:
#                 num_infecteds_average[sigma][val] = infecteds[k]
            
#             try:
#                 t_step_counts[sigma][val] += 1
#             except:
#                 t_step_counts[sigma][val] = 1
#         for k, val in enumerate(evolved_rhos):
#             if val > 1:
#                 print(seed, " ", k)
#             temp_rhos[k] = temp_rhos[k] + val/repeats

#         ax[0][0].plot(t_steps, evolved_rhos, label = fr"$\sigma$ = {sigma}", c = colours[i])
#         ax[0][1].plot(t_steps, infecteds, label = fr"$\sigma$ = {sigma}", c = colours[i])
    t_steps = list(average_rhos[sigma].keys())
#     print(t_steps)
    t_step_vals = list(t_step_counts[sigma].values())
#     print(t_step_vals)
    
#     print(average_rhos[sigma])
#     print(num_infecteds_average[sigma])
    
    rho_values = [val/num_sims for val in list(average_rhos[sigma].values())]
    num_infected_values = [val/num_sims for val in list(num_infecteds_average[sigma].values())]
#     max_infected = max(num_infected_values)
#     max_infected_ind = num_infected_values.index(max_infected)
    
#     t_steps_cleaned = [val/t_steps[max_infected_ind] for val in t_steps[0:max_infected_ind]]
#     rho_values_cleaned = [val for val in rho_values[0:max_infected_ind]]
#     infected_values_cleaned = [val for val in num_infected_values[0:max_infected_ind]]
#     for k, val in enumerate(rho_values):
#         if val > 1:
#             print(t_steps[k], rho_values[k], t_step_vals[k])
#     print(t_steps)
#     print(rho_values)
#     print(num_infected_values)
#     flop
    ax[0][0].plot(recording_mesh, rho_values, label = fr"$\sigma$ = {sigma}", c = colours[i])
    ax[0][1].plot(recording_mesh, num_infected_values, label = fr"$\sigma$ = {sigma}", c = colours[i])
#     ax[0].plot([k + 1 for k in range(len(temp_rhos))], temp_rhos, label = fr"$\sigma$ = {sigma}", c = colours[i])
    
print("Doing B")
average_rhos = {}
num_infecteds_average = {}
t_step_counts = {}
for i, zeta in enumerate(zetas):
    average_rhos[zeta] = {}
    num_infecteds_average[zeta] = {}
    
    t_step_counts[zeta] = {}
    num_sims = 0
    for j in range(repeats):
        solB.gillespie_simulation(seed, N, starting_infecteds, beta_max, c1, c2, zeta, eta, delta, alpha, gamma, mut_chance, t_max, resolution, S_increment)

        df = pd.read_csv(f"../data/gillespie_simulations_B/data_set{seed}.csv")
        
        seed += 1
        
        df = df[df["t_step"] <= t_max]
        
        evolved_rhos = [val for val in df["Average_rho"].values]
        t_steps = [val for val in df["t_step"].values]
        infecteds = [val for val in df["num_infected"].values]
        
        max_infecteds = max(infecteds)
        if max_infecteds > 1000:
            num_sims += 1
            inf_ind = infecteds.index(max_infecteds)
            t_peak = t_steps[inf_ind]
            recast_ts = [t/t_peak for t in t_steps[:inf_ind+1]]
    #         print(recast_ts[-1], t_peak)
            recast_infecteds = []
            recast_rhos = []
            rec_ind = 0
            for k, val in enumerate(recast_ts):
                if val >= recording_mesh[rec_ind]:
                    try:
                        average_rhos[zeta][recording_mesh[rec_ind]] += evolved_rhos[k]
                    except KeyError:
                        average_rhos[zeta][recording_mesh[rec_ind]] = evolved_rhos[k]

                    try:
                        num_infecteds_average[zeta][recording_mesh[rec_ind]] += infecteds[k]
                    except KeyError:
                        num_infecteds_average[zeta][recording_mesh[rec_ind]] = infecteds[k]

                    try:
                        t_step_counts[zeta][recording_mesh[rec_ind]] += 1
                    except KeyError:
                        t_step_counts[zeta][recording_mesh[rec_ind]] = 1


                    rec_ind += 1
#         df = df[df["t_step"] <= t_max]
        
#         seed += 1
        
#         evolved_rhos = df["Average_rho"].values
#         t_steps = df["t_step"].values
#         infecteds = df["num_infected"].values
#         for k, val in enumerate(evolved_rhos):
#             if val > 1:
#                 print(seed, t_steps[k], val, infecteds[k])
#         for k, val in enumerate(t_steps):
#             try:
#                 average_rhos[zeta][val] += evolved_rhos[k]
#             except:
#                 average_rhos[zeta][val] = evolved_rhos[k]
            
#             try:
#                 num_infecteds_average[zeta][val] += infecteds[k]
#             except:
#                 num_infecteds_average[zeta][val] = infecteds[k]
            
#             try:
#                 t_step_counts[zeta][val] += 1
#             except:
#                 t_step_counts[zeta][val] = 1
                
    t_steps = list(average_rhos[zeta].keys())
#     print(t_steps)
    t_step_vals = list(t_step_counts[zeta].values())
#     print(t_step_vals)
    
#     print(average_rhos[sigma])
#     print(num_infecteds_average[sigma])
    
    rho_values = [val/num_sims for val in list(average_rhos[zeta].values())]
    num_infected_values = [val/num_sims for val in list(num_infecteds_average[zeta].values())]
    
#     t_steps = list(average_rhos[zeta].keys())
#     t_step_vals = list(t_step_counts[zeta].values())
#     rho_values = [val/t_step_vals[k] for k, val in enumerate(list(average_rhos[zeta].values()))]
#     num_infected_values = [val/t_step_vals[k] for k, val in enumerate(list(num_infecteds_average[zeta].values()))]
#     max_infected = max(num_infected_values)
#     max_infected_ind = num_infected_values.index(max_infected)
    
#     t_steps_cleaned = [val/t_steps[max_infected_ind] for val in t_steps[0:max_infected_ind]]
#     rho_values_cleaned = [val for val in rho_values[0:max_infected_ind]]
#     infected_values_cleaned = [val for val in num_infected_values[0:max_infected_ind]]
#     for k, val in enumerate(rho_values):
#         if val > 1:
#             print(t_steps[k], rho_values[k], t_step_vals[k])
        
    ax[1][0].plot(recording_mesh, rho_values, label = fr"$\zeta$ = {zeta}", c = colours[i])
    ax[1][1].plot(recording_mesh, num_infected_values, label = fr"$\zeta$ = {zeta}", c = colours[i])
#     ax[1][0].plot(t_steps_cleaned, rho_values_cleaned, label = fr"$\zeta$ = {zeta}", c = colours[i])
#     ax[1][1].plot(t_steps_cleaned, infected_values_cleaned, label = fr"$\zeta$ = {zeta}", c = colours[i])
#         for k, val in enumerate(evolved_rhos):
#             if val > 1:
#                 print(seed, " ", k)
#             temp_rhos[k] = temp_rhos[k] + val/repeats
#     ax[1].plot([k + 1 for k in range(len(temp_rhos))], temp_rhos, label = fr"$\zeta$ = {zeta}", c = colours[i])

# ax[0].set_ylim([0,1])
# ax[1].set_ylim([0,1])

print("Doing C")
average_rhos = {}
num_infecteds_average = {}
t_step_counts = {}
for i, zeta in enumerate(zetas):
    average_rhos[zeta] = {}
    num_infecteds_average[zeta] = {}
    
    t_step_counts[zeta] = {}
    num_sims = 0
    for j in range(repeats):
        solC.gillespie_simulation(seed, N, starting_infecteds, beta_max, c1, c2, zeta, eta, delta, alpha, gamma, mut_chance, t_max, resolution, S_increment)

        df = pd.read_csv(f"../data/gillespie_simulations_C/data_set{seed}.csv")
        
        seed += 1
        
        df = df[df["t_step"] <= t_max]
        
        evolved_rhos = [val for val in df["Average_rho"].values]
        t_steps = [val for val in df["t_step"].values]
        infecteds = [val for val in df["num_infected"].values]
        
        max_infecteds = max(infecteds)
        if max_infecteds > 1000:
            num_sims += 1
            inf_ind = infecteds.index(max_infecteds)
            t_peak = t_steps[inf_ind]
            recast_ts = [t/t_peak for t in t_steps[:inf_ind+1]]
    #         print(recast_ts[-1], t_peak)
            recast_infecteds = []
            recast_rhos = []
            rec_ind = 0
            for k, val in enumerate(recast_ts):
                if val >= recording_mesh[rec_ind]:
                    try:
                        average_rhos[zeta][recording_mesh[rec_ind]] += evolved_rhos[k]
                    except KeyError:
                        average_rhos[zeta][recording_mesh[rec_ind]] = evolved_rhos[k]

                    try:
                        num_infecteds_average[zeta][recording_mesh[rec_ind]] += infecteds[k]
                    except KeyError:
                        num_infecteds_average[zeta][recording_mesh[rec_ind]] = infecteds[k]


                    rec_ind += 1
#         df = df[df["t_step"] <= t_max]
        
#         seed += 1
        
#         evolved_rhos = df["Average_rho"].values
#         t_steps = df["t_step"].values
#         infecteds = df["num_infected"].values
#         for k, val in enumerate(evolved_rhos):
#             if val > 1:
#                 print(seed, t_steps[k], val, infecteds[k])
#         for k, val in enumerate(t_steps):
#             try:
#                 average_rhos[zeta][val] += evolved_rhos[k]
#             except:
#                 average_rhos[zeta][val] = evolved_rhos[k]
            
#             try:
#                 num_infecteds_average[zeta][val] += infecteds[k]
#             except:
#                 num_infecteds_average[zeta][val] = infecteds[k]
            
#             try:
#                 t_step_counts[zeta][val] += 1
#             except:
#                 t_step_counts[zeta][val] = 1
                
    t_steps = list(average_rhos[zeta].keys())
#     print(t_steps)
    t_step_vals = list(t_step_counts[zeta].values())
#     print(t_step_vals)
    
#     print(average_rhos[sigma])
#     print(num_infecteds_average[sigma])
    
    rho_values = [val/num_sims for val in list(average_rhos[zeta].values())]
    num_infected_values = [val/num_sims for val in list(num_infecteds_average[zeta].values())]
#     print(colours[i])
    ax[2][0].plot(recording_mesh, rho_values, label = fr"$\zeta$ = {zeta}", c = colours[i])
    ax[2][1].plot(recording_mesh, num_infected_values, label = fr"$\zeta$ = {zeta}", c = colours[i])
    
ax[0][0].set_title("Model A")
ax[1][0].set_title("Model B")
ax[2][0].set_title("Model C")
ax[0][0].legend(loc = "upper left")
ax[1][0].legend(loc = "upper left")
ax[2][0].legend(loc = "upper left")

ax[0][0].set_ylim([0,0.4])
ax[1][0].set_ylim([0,0.4])
ax[2][0].set_ylim([0,0.4])

ax[0][1].set_ylim([0,70000])
ax[1][1].set_ylim([0,70000])
ax[2][1].set_ylim([0,70000])

ax[0][0].text(0.03, 1.010,labels[0][0], transform=ax[0][0].transAxes)
ax[1][0].text(0.03, 1.010,labels[1][0], transform=ax[1][0].transAxes)
ax[0][1].text(0.03, 1.010,labels[0][1], transform=ax[0][1].transAxes)
ax[1][1].text(0.03, 1.010,labels[1][1], transform=ax[1][1].transAxes)
ax[2][0].text(0.03, 1.010,labels[2][0], transform=ax[2][0].transAxes)
ax[2][1].text(0.03, 1.010,labels[2][1], transform=ax[2][1].transAxes)

ax[0][0].set_ylabel(r"Average detection avoidance, $\rho$")
ax[0][1].set_ylabel("Average number of infected individuals")
ax[1][1].set_xlabel(r"Time to peak infection, $\frac{t}{t_{peak}}$")
ax[0][1].set_xlabel(r"Time to peak infection, $\frac{t}{t_{peak}}$")
ax[2][1].set_xlabel(r"Time to peak infection, $\frac{t}{t_{peak}}$")

plt.savefig("../outputs/Gillespie_Models_imperfect_large_pop.png", bbox_inches = "tight")
plt.savefig("../outputs/Gillespie_Models_imperfect_large_pop.pdf", bbox_inches = "tight")
plt.close()