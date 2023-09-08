# distutils: language = c++
cimport cython
import numpy as np
from solvers cimport Quick_solver

cdef class PySolver:
    
    cdef Quick_solver cpp_solver
    
    def __cinit__(self):
        self.cpp_solver = Quick_solver()
        return
    
    
    def gillespie_simulation(self, int seed, int N, int starting_infecteds, float beta_max, float c1, float c2, float zeta, float eta, float delta, float alpha, float gamma, float mut_chance, float t_max, int resolution, int S_increment):
        
        self.cpp_solver.run_gillespie_simulation(seed, N, starting_infecteds, beta_max, c1, c2, zeta, eta, delta, alpha, gamma, mut_chance, t_max, resolution, S_increment)
        
        return 

