cdef extern from "wrapper.cpp":
    pass

cdef extern from "wrapper.h" namespace "solvers":
    cdef cppclass Quick_solver:
        Quick_solver() except +
        
        void run_gillespie_simulation(int, int, int, float, float, float, float, float, float, float, float, float, float, int, int)
                