#ifndef WRAPPER_H
#define WRAPPER_H

    namespace solvers{

        class Quick_solver{

            public:
                Quick_solver();
                
                ~Quick_solver();
                
                float fastmax(float a, float b);
                float fastmin(float a, float b);

                void create_proportions(float* output_array, float* input_array, int resolution);
                void beta_func(float* betas, float *rhos, float beta_max, float c1, float c2, int resolution);
                void infection_rate_function(float* infection_rate_array, int S, int* Qs, int* As, int* Us, float* betas, float delta, int resolution);
                void recovery_rate_function(float* recovery_rate_array, int* Qs, int* As, int* Us, float gamma, int resolution);
                void mortality_rate_function(float* mortality_rate_array, int* Qs, int* As, int* Us, float alpha, int resolution);
                void run_gillespie_simulation(int seed, int N, int starting_infecteds, float beta_max, float c1, float c2, float sigma, float eta, float delta, float alpha, float gamma, float mut_chance, float t_max, int resolution, int S_increment);
                
        };
    }


#endif