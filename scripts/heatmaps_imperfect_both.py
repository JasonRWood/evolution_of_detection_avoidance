import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
import math

sys.path.append("../src/model_A")

import runnerA
import fitness_functions_model_A as ffmA

sys.path.append("../src/model_B")

import runnerB
import fitness_functions_model_B as ffmB

sys.path.append("../src/model_C")

import runnerC
import fitness_functions_model_C as ffmC

solA = runnerA.PySolver()
solB = runnerB.PySolver()
solC = runnerC.PySolver()


b = 2.0
q = 0.1
d = 0.1
alpha = 0.5
gamma = 0.5

delta = 0.1
eta = 0.9

labels = ["(A)", "(B)", "(C)"]

c1 = 0.20

c2 = 4.0
seed = 31
# rho_init = 80

# sigma = 0.5

sigmas = [1.0*i/100 for i in range(101)]

# delta = 0.5

resolution = 150
S_density = I_F_density = 3.0
I_Q_density = 0.0
trait_space_length = 101

evo_steps = 2000

rho_max = 1.0
beta_max = 1.0

depth_max = 3

scenario = 3

#Specific Parameters 
sigma = 0.75 #For model A
zeta = 3.0 #For model B

output_mat = []

c1_max = 1.0
c1_min = 0.0

c2_max = 10
c2_min = -10

c1s = [c1_min + (c1_max - c1_min)*i/resolution for i in range(resolution + 1)]

c2s = [c2_min + (c2_max - c2_min)*i/resolution for i in range(resolution + 1)]

intermediate_vals_c1s_A = []
maximisers_c1s_A = []
minimisers_c1s_A = []
repeller_vals_c1s_A = []

intermediate_vals_c2s_A = []
maximisers_c2s_A = []
minimisers_c2s_A = []
repeller_vals_c2s_A = []

for c1 in c1s:
    for c2 in c2s:
        attractors, repellers = ffmA.find_sing_strats(solA, 0, depth_max, 10, 0, 1, beta_max, c1, c2, b, q, d, alpha, gamma, delta, eta, sigma, seed, scenario)
        
        if len(repellers) > 0:
            repeller_vals_c1s_A.append(c1)
            repeller_vals_c2s_A.append(c2)
            
        else:
            if attractors[0] == 0.0:
                minimisers_c1s_A.append(c1)
                minimisers_c2s_A.append(c2)
            elif attractors[0] == 1.0:
                maximisers_c1s_A.append(c1)
                maximisers_c2s_A.append(c2)
            else:
                intermediate_vals_c1s_A.append(c1)
                intermediate_vals_c2s_A.append(c2)
                    
intermediate_vals_c1s_B = []
maximisers_c1s_B = []
minimisers_c1s_B = []
repeller_vals_c1s_B = []

intermediate_vals_c2s_B = []
maximisers_c2s_B = []
minimisers_c2s_B = []
repeller_vals_c2s_B = []

for c1 in c1s:
    temp_row = []
    for c2 in c2s:
        attractors, repellers = ffmB.find_sing_strats(solB, 0, depth_max, 10, 0, 1, beta_max, c1, c2, b, q, d, alpha, gamma, delta, eta, zeta, seed, scenario)
        
        if len(repellers) > 0:
            repeller_vals_c1s_B.append(c1)
            repeller_vals_c2s_B.append(c2)
            
        else:
            if attractors[0] == 0.0:
                minimisers_c1s_B.append(c1)
                minimisers_c2s_B.append(c2)
            elif attractors[0] == 1.0:
                maximisers_c1s_B.append(c1)
                maximisers_c2s_B.append(c2)
            else:
                intermediate_vals_c1s_B.append(c1)
                intermediate_vals_c2s_B.append(c2)
                
intermediate_vals_c1s_C = []
maximisers_c1s_C = []
minimisers_c1s_C = []
repeller_vals_c1s_C = []

intermediate_vals_c2s_C = []
maximisers_c2s_C = []
minimisers_c2s_C = []
repeller_vals_c2s_C = []

for c1 in c1s:
    temp_row = []
    for c2 in c2s:
        attractors, repellers = ffmC.find_sing_strats(solC, 0, depth_max, 10, 0, 1, beta_max, c1, c2, b, q, d, alpha, gamma, delta, eta, zeta, seed, scenario)
        
        if len(repellers) > 0:
            repeller_vals_c1s_C.append(c1)
            repeller_vals_c2s_C.append(c2)
            
        else:
            if attractors[0] == 0.0:
                minimisers_c1s_C.append(c1)
                minimisers_c2s_C.append(c2)
            elif attractors[0] == 1.0:
                maximisers_c1s_C.append(c1)
                maximisers_c2s_C.append(c2)
            else:
                intermediate_vals_c1s_C.append(c1)
                intermediate_vals_c2s_C.append(c2)
    
fig = plt.figure(figsize = (15, 5))
gs = fig.add_gridspec(1,3)
ax = [fig.add_subplot(gs[i]) for i in range(3)]

ax[0].scatter(intermediate_vals_c2s_A, intermediate_vals_c1s_A, s = 1.25, label = "Intermediate value")
ax[0].scatter(minimisers_c2s_A, minimisers_c1s_A, s = 1.25, label = "Minimiser")
ax[0].scatter(maximisers_c2s_A, maximisers_c1s_A, s = 1.25, label = "Maximiser")
ax[0].scatter(repeller_vals_c2s_A, repeller_vals_c1s_A, s = 1.25, label = "Repeller")

ax[0].set_xlim([c2s[0], c2s[-1]])
ax[0].set_ylim([c1s[0], c1s[-1]])
ax[0].set_title("Model A")
ax[0].set_xlabel(r"Trade-off Curvature, $c_2$")
ax[0].set_ylabel(r"Trade-off Cost, $c_1$")
ax[0].text(0.03, 1.02,labels[0], transform=ax[0].transAxes)

ax[1].scatter(intermediate_vals_c2s_B, intermediate_vals_c1s_B, s = 1.25, label = "Intermediate value")
ax[1].scatter(minimisers_c2s_B, minimisers_c1s_B, s = 1.25, label = "Minimiser")
ax[1].scatter(maximisers_c2s_B, maximisers_c1s_B, s = 1.25, label = "Maximiser")
ax[1].scatter(repeller_vals_c2s_B, repeller_vals_c1s_B, s = 1.25, label = "Repeller")

ax[1].set_xlim([c2s[0], c2s[-1]])
ax[1].set_ylim([c1s[0], c1s[-1]])
ax[1].set_title("Model B")
ax[1].set_xlabel(r"Trade-off Curvature, $c_2$")
ax[1].set_ylabel(r"Trade-off Cost, $c_1$")
ax[1].text(0.03, 1.02,labels[1], transform=ax[1].transAxes)

ax[0].legend(loc = "upper left")

ax[2].scatter(intermediate_vals_c2s_C, intermediate_vals_c1s_C, s = 1.25, label = "Intermediate value")
ax[2].scatter(minimisers_c2s_C, minimisers_c1s_C, s = 1.25, label = "Minimiser")
ax[2].scatter(maximisers_c2s_C, maximisers_c1s_C, s = 1.25, label = "Maximiser")
ax[2].scatter(repeller_vals_c2s_C, repeller_vals_c1s_C, s = 1.25, label = "Repeller")

ax[2].set_xlim([c2s[0], c2s[-1]])
ax[2].set_ylim([c1s[0], c1s[-1]])
ax[2].set_title("Model C")
ax[2].set_xlabel(r"Trade-off Curvature, $c_2$")
ax[2].set_ylabel(r"Trade-off Cost, $c_1$")
ax[2].text(0.03, 1.02,labels[2], transform=ax[2].transAxes)

plt.savefig(f"../outputs/heatmap_imprerfect_both_with_C.png", bbox_inches = "tight")
plt.savefig(f"../outputs/heatmap_imprerfect_both_with_C.pdf", bbox_inches = "tight")
plt.close()