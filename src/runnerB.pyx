# distutils: language = c++
from solvers cimport Quick_solver

cdef class PySolver:
    
    cdef Quick_solver cpp_solver
    
    def __cinit__(self):
        self.cpp_solver = Quick_solver()
        return
    
    def ad_dyn(self, float beta_max, float rho_max, float b, float q, float d, float alpha_Q, float gamma_Q, float alpha_F, float gamma_F, float c1, float c2, float beta_Q, float beta_F, float delta_Q, float delta_F, float eta, float lam, float epsilon, float delta, int seed, int rho_init, float S_density = 4.0, float I_Q_1_density = 1.0, float I_Q_2_density = 1.0, float I_F_1_density = 1.0, float I_F_2_density = 1.0, int evo_steps = 1000, int trait_space_length = 101):
        
        self.cpp_solver.ad_dyn(beta_max, rho_max, b, q, d, alpha_Q, gamma_Q, alpha_F, gamma_F, c1, c2, beta_Q, beta_F, delta_Q, delta_F, eta, lam, epsilon, delta, seed, rho_init, S_density, I_Q_1_density, I_Q_2_density, I_F_1_density, I_F_2_density, evo_steps, trait_space_length)
        
        return
    
    def eco_dynamics(self, float beta_value, float rho_value, float b, float q, float d, float alpha_Q, float gamma_Q, float alpha_F, float gamma_F, float c1, float c2, float beta_Q, float beta_F, float delta_Q, float delta_F, float eta, float lam, float epsilon, float delta, int seed, float S_density = 4.0, float I_Q_1_density = 1.0, float I_Q_2_density = 1.0, float I_F_1_density = 1.0, float I_F_2_density = 1.0):

        cdef float s, i_q_1, i_q_2, i_f_1, i_f_2, r
        
        self.cpp_solver.eco_dynamics(&s, &i_q_1, &i_q_2, &i_f_1, &i_f_2, &r, beta_value, rho_value, b, q, d, alpha_Q, gamma_Q, alpha_F, gamma_F, c1, c2, beta_Q, beta_F, delta_Q, delta_F, eta, lam, epsilon, delta, seed, S_density, I_Q_1_density, I_Q_2_density, I_F_1_density, I_F_2_density)

        y = [s,i_q_1,i_q_2,i_f_1,i_f_2]
        
        return y

