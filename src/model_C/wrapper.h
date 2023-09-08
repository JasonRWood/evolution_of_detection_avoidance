    #ifndef WRAPPER_H
    #define WRAPPER_H

    namespace solvers{

        class Quick_solver{

            public:
                float beta_max, alpha_max, sigma_max, b, q, d, rho, eta, gamma, lambda, c1, c2, hyper;
                int seed;
                Quick_solver();
                Quick_solver(float beta_max, float alpha_max, float sigma_max, float b, float q, float d, float rho, float eta, float gamma, float lambda, float c1, float c2, float hyper, int seed);
                ~Quick_solver();
                
                float fastmax(float a, float b);
                float fastmin(float a, float b);
                int dynamics(float* dydt, float* y, int num_parasites, int* rho_inds,  float b, float q, float d, float *beta, float alpha, float gamma, float* rho, float delta_Q, float eta, float delta);
                int rkck(float *y, float *dydt, float *yout, float *yerr, float h, int num_parasites, int* rho_inds,  float b, float q, float d, float* beta, float alpha, float gamma, float* rho, float delta_Q, float eta, float delta);
                int rkqs(float *y, float *dydt, float *h,float *hnext,float *yscale, int num_parasites, int* rho_inds,  float b, float q, float d, float* beta, float alpha, float gamma, float* rho, float delta_Q, float eta, float delta, float* t);
                int ad_dyn(float beta_max, float rho_max, float b, float q, float d, float alpha, float gamma, float c1, float c2, float delta_Q, float eta, float delta, int seed, int rho_init, float S_density, float Q_density, float Adensity, float U_density, float I_density, int evo_steps, int trait_space_length);
                void eco_dynamics(float* S, float* Q, float* A, float* U, float* I, float* R, float beta_value, float rho_value, float b, float q, float d, float alpha, float gamma, float c1, float c2, float delta_Q, float eta, float delta, int seed, float S_density, float Q_density, float A_density, float U_density, float I_density);
        };
    }


    #endif