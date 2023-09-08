cdef extern from "wrapper.cpp":
    pass

cdef extern from "wrapper.h" namespace "solvers":
    cdef cppclass Quick_solver:
        Quick_solver() except +
        Quick_solver(float, float, float, float, float, float, float, float, float, float, float, float, float, int) except +
        int ad_dyn(float, float, float, float, float, float, float, float, float, float, float, float, float, float, float, float, float, float, float, int, int, float, float, float, float, float, int, int)
        int ad_dyn(float , float , float , float , float , float , float , float , float , float , float , float , int , int , float , float , float , float , int , int )
        void eco_dynamics(float*, float*, float*, float*, float*, float, float, float, float, float, float, float, float, float, float, float, float, int, float, float, float, float)