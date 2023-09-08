# distutils: language = c++
from solvers cimport Quick_solver

cdef class PySolver:
    
    cdef Quick_solver cpp_solver
    
    def __cinit__(self):
        self.cpp_solver = Quick_solver()
        return
    
    def ad_dyn(self, float beta_max, float rho_max, float b, float q, float d, float alpha, float gamma, float c1, float c2, float delta_Q, float eta, float sigma, int seed, int rho_init, float S_density = 4.0, float Q_density = 1.0, float A_density = 1.0, float U_density = 1.0, int evo_steps = 1000, int trait_space_length = 101):
        
        self.cpp_solver.ad_dyn(beta_max, rho_max, b, q, d, alpha, gamma, c1, c2, delta_Q, eta, sigma, seed, rho_init, S_density, Q_density, A_density, U_density, evo_steps, trait_space_length)
        
        return
    
    def eco_dynamics(self, float beta_value, float rho_value, float b, float q, float d, float alpha, float gamma, float c1, float c2, float delta_Q, float eta, float sigma, int seed, float S_density = 1.0, float Q_density = 1.0, float A_density = 1.0, float U_density = 1.0):

        cdef float S, Q, A, U, R
        
        self.cpp_solver.eco_dynamics(&S, &Q, &A, &U, &R, beta_value, rho_value, b, q, d, alpha, gamma, c1, c2, delta_Q, eta, sigma, seed, S_density, Q_density, A_density, U_density)
        
        y = [S, Q, A, U, R]
        
        return y

