#include <iostream>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <time.h>
#include <cmath>
#include <unistd.h>
#include <iomanip> 
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <wrapper.h>
using namespace std;


namespace solvers{

    Quick_solver::Quick_solver() {};
    
    Quick_solver::~Quick_solver() {};

    float Quick_solver::fastmax(float a, float b){
        if (a > b){
            return a;
        }
        else{
            return b;
        }
    }

    float Quick_solver::fastmin(float a, float b){
        if (a < b){
            return a;
        }
        else{
            return b;
        }
    }

    void Quick_solver::create_proportions(float* output_array, float* input_array, int resolution){

        float temp_sum = 0.0;

        for (int i = 0; i < resolution + 1; i++){

            temp_sum += input_array[i];
        }

        output_array[0] = input_array[0]/temp_sum;
        for (int i = 1; i < resolution + 1; i++){
            output_array[i] = output_array[i - 1] + input_array[i]/temp_sum;
        }
        return;
    }

    void Quick_solver::beta_func(float* betas, float *rhos, float beta_max, float c1, float c2, int resolution){
        
        float rho_max = 1.0;
        for(int i=0;i<resolution + 1;i++){
            if (c2 == 0.0){
                betas[i] = beta_max*((1-c1) + c1*((rho_max - rhos[i])/rho_max));
            }
            else{
                betas[i] = beta_max*((1-c1) + c1*((1 - exp((rho_max - rhos[i])*c2/rho_max))/(1 - exp(c2))));
            }
        }
        return;
    }

    void Quick_solver::infection_rate_function(float* infection_rate_array, int S, int* Qs, int* As, int* Us, float* betas, float delta, int resolution){

        for (int i = 0; i< resolution + 1; i++){
            float Infection_rate;
            Infection_rate = betas[i]*(delta*Qs[i] + As[i] + Us[i])*S;
            infection_rate_array[i] = Infection_rate;

        }
        return;
    }

    void Quick_solver::recovery_rate_function(float* recovery_rate_array, int* Qs, int* As, int* Us, float gamma, int resolution){

        for (int i = 0; i< resolution + 1; i++){
            float recovery_rate;
            recovery_rate = gamma*(Qs[i] + As[i] + Us[i]);
            recovery_rate_array[i] = recovery_rate;

        }
        return;
    }

    void Quick_solver::mortality_rate_function(float* mortality_rate_array, int* Qs, int* As, int* Us, float alpha, int resolution){

        for (int i = 0; i< resolution + 1; i++){
            float mortality_rate;
            mortality_rate = alpha*(Qs[i] + As[i] + Us[i]);
            mortality_rate_array[i] = mortality_rate;

        }
        return;
    }

    void Quick_solver::run_gillespie_simulation(int seed, int N, int starting_infecteds, float beta_max, float c1, float c2, float sigma, float eta, float delta, float alpha, float gamma, float mut_chance, float t_max, int resolution, int S_increment){
        
        srand(seed);
        string seed_str;
        seed_str = std::to_string(seed);
        ofstream output_file;

        output_file.open("../data/gillespie_simulations_A/data_set" + seed_str +  ".csv");
        output_file << "t_step,Average_rho,num_infected\n";

        float* rhos = new float[resolution + 1];
        float* betas = new float[resolution + 1];
        float rho_max = 1.0;
        int i = 0;

        int* proc_ind = new int[1];
        float* average_rho = new float[1];

        // int fail_ind = -1;
        // float* rates_vectors = new float[max_steps];
        // float* evolved_rhos = new float[max_steps]; 

        for (i = 0; i < resolution + 1; i++){
            rhos[i] = (rho_max*(i))/(resolution);
        }

        beta_func(betas, rhos, beta_max, c1, c2, resolution);

        // for(i = 0; i < resolution + 1; i++){
        //     std::cout << "For index " << i << " the beta value is " << betas[i] << "\n";
        // }
        // if(beta_max > 0){
        //     return;
        // }
        int S = N - starting_infecteds;
        int* Qs = new int[resolution + 1];
        int* As = new int[resolution + 1];
        int* Us = new int[resolution + 1];

        float* pops = new float[3];
        float* props = new float[3];

        for(i = 0; i < resolution + 1; i++){
            Qs[i] = 0;
            As[i] = 0;
            Us[i] = 0;
        }
        Us[0] = starting_infecteds;

        float t_curr = 0;
        float t_check;
        float t_check_increment;

        if (N > 10000){
            t_check = 0.00001;
            t_check_increment = 0.00001;
        }
        else{
            t_check = 0.01;
            t_check_increment = 0.01;
        }
        

        average_rho[0] = 0;
        while(t_curr < t_max){
            
            // output_file << k << ",";

            float* infection_rate_array = new float[resolution + 1];
            float* mortality_rate_array = new float[resolution + 1];
            float* recovery_rate_array = new float[resolution + 1];

            infection_rate_function(infection_rate_array, S, Qs, As, Us, betas, delta, resolution);
            recovery_rate_function(recovery_rate_array, Qs, As, Us, gamma, resolution);
            mortality_rate_function(mortality_rate_array, Qs, As, Us, alpha, resolution);

            float inf_sum = 0, rec_sum = 0, mort_sum = 0;

            for (i = 0; i < resolution + 1; i++){
                inf_sum += infection_rate_array[i];
                rec_sum += recovery_rate_array[i];
                mort_sum += mortality_rate_array[i];
            }

            float* total_rate_array = new float[3];
            total_rate_array[0] = inf_sum;
            total_rate_array[1] = rec_sum;
            total_rate_array[2] = mort_sum;

            float rand_stored_proc;
            rand_stored_proc = (float(rand())/float((RAND_MAX)));

            proc_ind[0] = 0;

            float* total_prop_array = new float[3];

            total_prop_array[0] = inf_sum/(inf_sum + rec_sum + mort_sum);
            total_prop_array[1] = (inf_sum + rec_sum)/(inf_sum + rec_sum + mort_sum);
            total_prop_array[2] = 1.0;
            // create_proportions(total_prop_array, total_rate_array, 2);

            float rate = (inf_sum + rec_sum + mort_sum);
            float rand_stored_time;
            rand_stored_time = (float(rand())/float((RAND_MAX)));
            float t_inc = 1/rate*log(1/rand_stored_time);
            
                
            if (t_curr + t_inc > t_check){
                //Record
                t_curr = t_check;
                int num_infected = 0;
                for (i = 0; i < resolution + 1; i++){
                    num_infected += Qs[i];
                    num_infected += As[i];
                    num_infected += Us[i];
                }
                output_file << t_check << "," << average_rho[0] << "," << num_infected << "\n";
                delete recovery_rate_array;
                delete infection_rate_array;
                delete mortality_rate_array;
                t_check += t_check_increment;
            }
            else{
                t_curr += t_inc;
                while(total_prop_array[proc_ind[0]] < rand_stored_proc){
                    proc_ind[0] += 1;
                }

                // if (proc_ind[0] > 3){
                //     std::cout << "proc_ind is " << proc_ind[0] << "\n";
                //     for(i = 0; i < 3; i++){
                //         std::cout << "total_rate_array[" << i << "] = " << total_rate_array[i] << "\n";
                //         std::cout << "total_prop_array[" << i << "] = " << total_prop_array[i] << "\n";
                //     }
                // }
                delete total_prop_array;
                delete total_rate_array;

                if(proc_ind[0] == 0){
                        // Do infection
                        S += -1;
                        float rand_stored_infector = (float(rand())/float((RAND_MAX)));
                        float* infection_proportion_array = new float[resolution + 1];
                        create_proportions(infection_proportion_array, infection_rate_array, resolution);

                        int infector_ind = 0;
                        while(infection_proportion_array[infector_ind] < rand_stored_infector){
                            infector_ind += 1;
                        }

                        delete infection_proportion_array;

                        float chosen_propensity = infection_rate_array[infector_ind];

                        float t_inc = 0;
                        // t_inc = 
                        // output_file << infection_rate_array[infector_ind] << ",";

                        if (float(rand())/float((RAND_MAX)) < mut_chance){
                            if (infector_ind != 0 && infector_ind != resolution){

                                float aware_rand = (float(rand())/float((RAND_MAX)));
                                float mut_rand = (float(rand())/float((RAND_MAX)));

                                if (sigma*rhos[infector_ind] + (1 - sigma) > aware_rand){
                                    if (mut_rand > 0.5){
                                        Us[infector_ind - 1] += 1;
                                    }
                                    else{
                                        Us[infector_ind + 1] += 1;
                                    }
                                }
                                else if (sigma*rhos[infector_ind] + (1 - sigma) + sigma*(1 - rhos[infector_ind])*eta > aware_rand){
                                    if (mut_rand > 0.5){
                                        Qs[infector_ind - 1] += 1;
                                    }
                                    else{
                                        Qs[infector_ind + 1] += 1;
                                    }
                                }
                                else{
                                    if (mut_rand > 0.5){
                                        As[infector_ind - 1] += 1;
                                    }
                                    else{
                                        As[infector_ind + 1] += 1;
                                    }
                                }                     
                                
                                }
                            else if (infector_ind == 0){
                                float aware_rand = (float(rand())/float((RAND_MAX)));
                                if (sigma*rhos[infector_ind] + (1 - sigma) > aware_rand){
                                    Us[infector_ind + 1] += 1;
                                }
                                else if (sigma*rhos[infector_ind] + (1 - sigma) + sigma*(1 - rhos[infector_ind])*eta > aware_rand){
                                    Qs[infector_ind + 1] += 1;
                                }
                                else{
                                    As[infector_ind + 1] += 1;
                                }
                            }
                            else if (infector_ind == resolution){
                                float aware_rand = (float(rand())/float((RAND_MAX)));
                                if (sigma*rhos[infector_ind] + (1 - sigma) > aware_rand){
                                    Us[infector_ind - 1] += 1;
                                }
                                else if (sigma*rhos[infector_ind] + (1 - sigma) + sigma*(1 - rhos[infector_ind])*eta > aware_rand){
                                    Qs[infector_ind - 1] += 1;
                                }
                                else{
                                    As[infector_ind - 1] += 1;
                                }

                            } 

                            }
                        else{
                            float aware_rand = (float(rand())/float((RAND_MAX)));
                            if (sigma*rhos[infector_ind] + (1 - sigma) > aware_rand){
                                Us[infector_ind] += 1;
                            }
                            else if (sigma*rhos[infector_ind] + (1 - sigma) + sigma*(1 - rhos[infector_ind])*eta > aware_rand){
                                Qs[infector_ind] += 1;
                            }
                            else{
                                As[infector_ind] += 1;
                            }
                        }
                    }
                    if (proc_ind[0] == 1){
                        // Do recovery
                        S += S_increment;
                        float rand_stored_recoverer = (float(rand())/float((RAND_MAX)));
                        float* recovery_proportion_array = new float[resolution + 1];
                        create_proportions(recovery_proportion_array, recovery_rate_array, resolution);

                        int rec_ind = 0;
                        while(recovery_proportion_array[rec_ind] < rand_stored_recoverer){
                            rec_ind += 1;
                        }
                        
                        delete recovery_proportion_array;

                        // output_file << recovery_rate_array[rec_ind] << ",";
                        int Q = Qs[rec_ind], A = As[rec_ind], U = Us[rec_ind];
                        // float pops[3];
                        // float props[3];
                        pops[0] = Q;
                        pops[1] = A;
                        pops[2] = U;

                        props[0] = float(Q)/(Q + A + U);
                        props[1] = float(Q + A)/(Q + A + U);
                        props[2] = 1.0;

                        int choice_ind = 0;
                        float choice_rand = (float(rand())/float((RAND_MAX)));
                        while(props[choice_ind] < choice_rand){
                            choice_ind += 1;
                        }

                        if (choice_ind == 0){
                            Qs[rec_ind] += -1;
                        }
                        if (choice_ind == 1){
                            As[rec_ind] += -1;
                        }
                        if (choice_ind == 2){
                            Us[rec_ind] += -1;
                        }
                    }

                    if (proc_ind[0] == 2){
                        // Do mortality
                        S += S_increment;
                        float rand_stored_dies = (float(rand())/float((RAND_MAX)));
                        float* mortality_proportion_array = new float[resolution + 1];
                        create_proportions(mortality_proportion_array, mortality_rate_array, resolution);

                        int die_ind = 0;
                        while(mortality_proportion_array[die_ind] < rand_stored_dies){
                            die_ind += 1;
                        }

                        delete mortality_proportion_array;

                        // output_file << mortality_rate_array[die_ind] << ",";

                        int Q = Qs[die_ind], A = As[die_ind], U = Us[die_ind];

                        pops[0] = Q;
                        pops[1] = A;
                        pops[2] = U;

                        props[0] = float(Q)/(Q + A + U);
                        props[1] = float(Q + A)/(Q + A + U);
                        props[2] = 1.0;

                        int choice_ind = 0;
                        float choice_rand = (float(rand())/float((RAND_MAX)));
                        while(props[choice_ind] < choice_rand){
                            choice_ind += 1;
                        }

                        if (choice_ind == 0){
                            Qs[die_ind] += -1;
                        }
                        if (choice_ind == 1){
                            As[die_ind] += -1;
                        }
                        if (choice_ind == 2){
                            Us[die_ind] += -1;
                        }
                    }
                delete recovery_rate_array;
                delete infection_rate_array;
                delete mortality_rate_array;

                float total_pop = 0;
                for (i = 0; i < resolution + 1; i ++){
                    total_pop += Qs[i] + As[i] + Us[i];
                }
                if (total_pop == 0){
                    std::cout << "Infection went extinct" << endl;
                    // output_file << 0 << "," << proc_ind[0] << "\n";
                    // fail_ind = k;
                    break;
                }
                average_rho[0] = 0;
                for (i = 0; i < resolution + 1; i ++){
                    float pop_proportion = float(Qs[i] + As[i] + Us[i])/total_pop;
                    average_rho[0] += rhos[i]*pop_proportion;
                }
                
                // if (average_rho[0] > 1){
                //     for (i = 0; i < resolution + 1; i ++){
                //         std::cout << "Qs[" << i << "] = " << Qs[i] << "\n";
                //         std::cout << "As[" << i << "] = " << As[i] << "\n";
                //         std::cout << "Us[" << i << "] = " << Us[i] << "\n";
                //         std::cout << "rhos[" << i << "] = " << rhos[i] << "\n";
                //         std::cout << "betas[" << i << "] = " << betas[i] << "\n";
                // }
                // }
                // output_file << average_rho[0] << "," << proc_ind[0] << "\n";
                
                // delete average_rho;
            }
            
            }
       

        // if (fail_ind == -1){
        //     for (int k = 0; k < max_steps; k ++){
        //         output_file << k << "," << rates_vectors[k] << "," << evolved_rhos[k] << "\n";
        //     }
        // }
        // else{
        //     for (int k = 0; k < fail_ind; k ++){
        //         output_file << k << "," << rates_vectors[k] << "," << evolved_rhos[k] << "\n";
        //     }
        // }
        output_file.close();
        // delete rates_vectors;
        // delete evolved_rhos;
        delete rhos;
        delete betas;

        delete Qs;
        delete As;
        delete Us;

        delete pops;
        delete props;

        delete proc_ind;
        delete average_rho;
        return;
        }
    }
