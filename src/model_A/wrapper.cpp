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

float TINY = 1e-3;
float EPS = 1e-6;

// int trait_space_length = 201;

namespace solvers{
    
    Quick_solver::Quick_solver() {};
    Quick_solver::Quick_solver(float beta_max, float alpha_max, float sigma_max, float b, float q, float d, float rho, float eta, float gamma, float lambda, float c1, float c2, float hyper, int seed){
        this->beta_max = beta_max;
        this->alpha_max = alpha_max;
        this->sigma_max = sigma_max;
        this->b = b;
        this->q = q;
        this->d = d;
        this->rho = rho;
        this->eta = eta;
        this->gamma = gamma;
        this->lambda = lambda;
        this->c1 = c1;
        this->c2 = c2;
        this->hyper = hyper;
        this->seed = seed;

    };
    
    
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

            // float fastmin(float a, float b){
            //     return (a<b)?a:b;
            // }

    int Quick_solver::dynamics(float* dydt, float* y, int num_parasites, int* rho_inds,  float b, float q, float d, float *beta, float alpha, float gamma, float* rho, float delta_Q, float eta, float sigma){
        
        int i;
        float Isum = 0.0, N = 0.0;
        float infection_force_para[num_parasites], para_infection_force = 0.0;
        
        float S, R;
        float Q[num_parasites], A[num_parasites];
        float U[num_parasites];
        
        float Sdot, Rdot;
        float Q_dot[num_parasites], A_dot[num_parasites];
        float U_dot[num_parasites];
            
        
        S = y[0];
        R = y[3*num_parasites + 1];
        for(i=0;i<num_parasites;i++){
            Q[i] = y[i + 1];
            A[i] = y[i + 1 + num_parasites];
            U[i] = y[i + 1 + 2*num_parasites];
            
            Isum +=  Q[i] + A[i] + U[i];
            
            infection_force_para[i] = beta[rho_inds[i]]*(delta_Q*Q[i] + A[i] + U[i]);
            para_infection_force += infection_force_para[i];
        }
        
        N = Isum + S + R;

        Sdot =  b*(1 - q*N)*N - (para_infection_force + d)*S;
        
        for(i = 0;i<num_parasites;i++){

            Q_dot[i] = sigma*(1 - rho[rho_inds[i]])*eta*infection_force_para[i]*S - (d + alpha + gamma)*Q[i];
            
            A_dot[i] = sigma*(1 - rho[rho_inds[i]])*(1 - eta)*infection_force_para[i]*S - (d + alpha + gamma)*A[i];

            U_dot[i] = (sigma*rho[rho_inds[i]] + (1 - sigma))*infection_force_para[i]*S - (d + alpha + gamma)*U[i];

        }

        Rdot = gamma*Isum  -  d*R;
        dydt[0] = Sdot;

        for(i = 0;i<num_parasites;i++){

            dydt[i + 1] = Q_dot[i];
            dydt[i + 1 + num_parasites] = A_dot[i];
            dydt[i + 1 + 2*num_parasites] = U_dot[i];
        }
        dydt[3*num_parasites + 1] = Rdot;
        
        return 0;
    }

    int Quick_solver::rkck(float *y, float *dydt, float *yout, float *yerr, float h, int num_parasites, int* rho_inds,  float b, float q, float d, float* beta, float alpha, float gamma, float* rho, float delta_Q, float eta, float sigma)
    {
        int i;
        // float y_1[2 + 3*num_parasites ], y_2[2 + 3*num_parasites ], y_3[2 + 3*num_parasites ], y_4[2 + 3*num_parasites ], y_5[2 + 3*num_parasites ];
        float* y_1 = new float[2 + 3*num_parasites];
        float* y_2 = new float[2 + 3*num_parasites];
        float* y_3 = new float[2 + 3*num_parasites];
        float* y_4 = new float[2 + 3*num_parasites];
        float* y_5 = new float[2 + 3*num_parasites];
        float* y_temp = new float[2 + 3*num_parasites];
        // float y_temp[2 + 3*num_parasites ];
        static float b21=0.2,b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42=-0.9,b43=1.2,
        b51=-11.0/54.0,b52=2.5,b53=-70.0/27.0,b54=35.0/27.0,b61=1631.0/55296,
        b62=175.0/512.0,b63=575.0/13824.0,b64=44275.0/110592,b65=253.0/4096.0,
        c1=37.0/378.0,c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,dc5=-277.00/14336;
        float dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,dc4=c4-13525.0/55296.0,
        dc6=c6-0.25;
        
        for(i=0;i<(2 + 3*num_parasites);i++){
            y_temp[i] = y[i] + b21*h*dydt[i];
        }
        
        dynamics(y_1, y_temp, num_parasites, rho_inds, b, q, d, beta, alpha, gamma, rho, delta_Q, eta, sigma);

        for(i=0;i<(2 + 3*num_parasites);i++){
            y_temp[i] = y[i]+h*(b31*dydt[i]+b32*y_1[i]);
        }
        
        dynamics(y_2, y_temp, num_parasites, rho_inds, b, q, d, beta, alpha, gamma, rho, delta_Q, eta, sigma);
        
        for(i=0;i<(2 + 3*num_parasites);i++){
            y_temp[i]= y[i]+h*(b41*dydt[i]+b42*y_1[i]+b43*y_2[i]);
        }

        dynamics(y_3, y_temp, num_parasites, rho_inds, b, q, d, beta, alpha, gamma, rho, delta_Q, eta, sigma);
        
        for(i=0;i<(2 + 3*num_parasites);i++){
            y_temp[i]= y[i]+h*(b51*dydt[i]+b52*y_2[i]+b53*y_2[i]+b54*y_3[i]);
        }

        dynamics(y_4, y_temp, num_parasites, rho_inds, b, q, d, beta, alpha, gamma, rho, delta_Q, eta, sigma);
        
        for(i=0;i<(2 + 3*num_parasites);i++){
            y_temp[i]= y[i]+h*(b61*dydt[i]+b62*y_1[i]+b63*y_2[i]+b64*y_3[i]+b65*y_4[i]);
        }

        dynamics(y_5, y_temp, num_parasites, rho_inds, b, q, d, beta, alpha, gamma, rho, delta_Q, eta, sigma);
        
        for(i=0;i<(2 + 3*num_parasites);i++){
            yout[i]= y[i]+h*(c1*dydt[i]+c3*y_2[i]+c4*y_3[i]+c6*y_5[i]);
            yerr[i]= h*(dc1*dydt[i]+dc3*y_2[i]+dc4*y_3[i]+dc5*y_4[i]+dc6*y_5[i]);
        }

        delete y_1;
        delete y_2;
        delete y_3;
        delete y_4;
        delete y_5;
        delete y_temp;
        return 0;
    }

    int Quick_solver::rkqs(float *y, float *dydt, float *h,float *hnext,float *yscale, int num_parasites, int* rho_inds,  float b, float q, float d, float* beta, float alpha, float gamma, float* rho, float delta_Q, float eta, float sigma, float* t)
    {
        
        float y_temp[2 + 3*num_parasites], yerr[2 + 3*num_parasites];
        float htemp,errmax;
        int i;//, counter_internal;
        
        htemp= *h;
        
        rkck(y, dydt, y_temp, yerr, *h, num_parasites, rho_inds, b, q, d, beta, alpha, gamma, rho, delta_Q, eta, sigma);
        *hnext= *h;
        
        for(;;)
        {
            rkck(y, dydt, y_temp, yerr, *h, num_parasites, rho_inds, b, q, d, beta, alpha, gamma, rho, delta_Q, eta, sigma);
            errmax= 0.0;
            
            for(i=0;i<(2 + 3*num_parasites);i++){
                if (errmax < abs(yerr[i]/yscale[i])){
                    errmax = abs(yerr[i]/yscale[i]);
                }
                
            }
            errmax/= 1e-4;
            if(errmax<=1.0){ 
                break;
            }
            
            htemp= 0.9*(*h)*pow(errmax,-0.25);
            
            *h= (*h>=0.0 ? fastmax(htemp,0.1*(*h)) : fastmin(htemp,0.1*(*h)));
        }
        
        if(errmax > 1.89E-4) {
            *hnext= 0.9*(*h)*pow(errmax,-0.2);
        } 
        else {
            *hnext= 5.0*(*h);
        }
        
        for (i = 0; i<(2 + 3*num_parasites); i++){
            y[i] = y_temp[i];
        }

        return 0;
    }

    int Quick_solver::ad_dyn(float beta_max, float rho_max, float b, float q, float d, float alpha, float gamma, float c1, float c2, float delta_Q, float eta, float sigma, int seed, int rho_init, float S_density, float Q_density, float A_density, float U_density, int evo_steps, int trait_space_length){
        
        srand(seed);
        
        int tmax = 1000, i;
        float beta[trait_space_length], rho[trait_space_length];
        float t[1], tcheck[1];

        int num_parasites = 1, rho_inds[trait_space_length], ind,num_parasites_2;
        int num_parasites_cleaned;

        float* y = new float[2+ 3*trait_space_length];
        float* y_check = new float[2+ 3*trait_space_length];
        float* dydt = new float[2+ 3*trait_space_length];
        float* yScale = new float[2+ 3*trait_space_length];
        float* y_max = new float[2+ 3*trait_space_length];
        float* y_min = new float[2+ 3*trait_space_length];
        float* y_max_next = new float[2+ 3*trait_space_length];
        float* y_min_next = new float[2+ 3*trait_space_length];

        float discrep_check;
        int rho_inds_cleaned[trait_space_length];
        
        float Q_temp[trait_space_length], U_temp[trait_space_length];
        float A_temp[trait_space_length];
        
        int rho_inds_2[trait_space_length];
        
        float* y_temp = new float[2+ 3*trait_space_length ];
        
        float* para_density = new float[trait_space_length];
        float* para_cum_density = new float[trait_space_length];
        float* para_props = new float[trait_space_length];

        float h[1], hnext[1];

        int num_poss_outputs = 40000;
        int para_output_counter = 0;

        int para_output_rho_inds[num_poss_outputs];
        float para_output_Q_densities[num_poss_outputs], para_output_A_densities[num_poss_outputs], para_output_rho_vals[num_poss_outputs];
        float para_output_U_densities[num_poss_outputs];
        float para_output_host_density[num_poss_outputs];
        int para_output_evosteps[num_poss_outputs];
        float para_output_R_densities[num_poss_outputs];
        ind = 0;
        ofstream para_outputs, tracker_file;
        ofstream hyper_outputs;
        ofstream joint_outputs;

        for(i=0;i<trait_space_length;i++){
            rho[i] = (rho_max*(i))/(trait_space_length - 1);
        }

        for(i=0;i<trait_space_length;i++){
            if (c2 == 0.0){
                beta[i] = beta_max*((1-c1) + c1*((rho_max - rho[i])/rho_max));
            }
            else{
                beta[i] = beta_max*((1-c1) + c1*((1 - exp((rho_max - rho[i])*c2/rho_max))/(1 - exp(c2))));
            }
        }

        ofstream beta_values,rho_values;
        beta_values.open("../data/beta_vals.csv");
        rho_values.open("../data/rho_vals.csv");

        for (i = 0; i<trait_space_length;i++){
            beta_values << i << ",";
            
            rho_values << i << ",";
            
        }
        beta_values << "\n";
        rho_values << "\n";
        

        

        for(i=0;i<trait_space_length;i++){
            rho_values << rho[i] << ",";
            beta_values << beta[i] << ",";
        }

        beta_values << "\n";
        rho_values << "\n";

        for (i = 0; i<trait_space_length;i++){
            rho_values << beta[i]/beta_max << ",";
            
        }
        rho_values << "\n";

        beta_values.close();
        rho_values.close();
        

        y[0] = S_density/(b/q);
        y[1] = Q_density/(b/q);
        y[2] = A_density/(b/q);
        y[3] =  U_density/(b/q);
        
        y[4] = 0.0;
        rho_inds[0] = rho_init;
        
        num_parasites = 1;
        

        // cout << y[0] << endl;
        // cout << y[1] << endl;
        // cout << y[2] << endl;

        for (int evo_counter = 0; evo_counter < evo_steps; evo_counter++){
            
            
            t[0] = 0.0;
            h[0] = 1e-1;
            hnext[0] = 1e-1;

            float y_check_sum = 0.0;
            
            tcheck[0] = 10.0;
            for(i=0;i<(2 + 3*num_parasites );i++){
                y_check[i] = y[i];
                y_max[i] = y[i];
                y_min[i] = y[i];
                y_max_next[i] = y[i];
                y_min_next[i] = y[i];
                y_check_sum += y[i];
                
            }
            int step_count = 0;
            int check_count = 0;

            // cout << "Before dynamics " << endl;
            // cout << y[0] << endl;
            // for(i = 0;i<num_parasites;i++){
            //     cout << "Parasite " << i <<" has I_Q_1 density "<< y[i + 1] <<" has I_Q_2 density "<< y[i + 1 + num_parasites] << endl;
            //     cout << " and I_F_1 density " << y[i + 1 +2*num_parasites] << " and I_F_2 density " << y[i + 1 + 3*num_parasites] << endl;
            // }
            // cout << "With " << y[3*num_parasites + 1] << " recovered hosts" << endl;

            while (t[0] <= tmax){

                
                dynamics(dydt, y, num_parasites, rho_inds, b, q, d, beta, alpha, gamma, rho, delta_Q, eta, sigma);

                /* Adjust the step size to maintain accuracy */
                for (i=0; i<(2 + 3*num_parasites ); i++){
                    yScale[i]=fabs(y[i])+fabs(dydt[i]*(*h))+EPS;
                }
                
                
                rkqs(y, dydt, h,hnext,yScale, num_parasites, rho_inds, b, q, d, beta, alpha, gamma, rho, delta_Q, eta, sigma, t);
                
                // cout << "During dynamics " << endl;
                // cout << "Host has density " <<  y[0] << " and" <<  endl;
                // for(i = 0;i<num_parasites;i++){
                //     cout << "Parasite " << i <<" has I_Q_1 density "<< y[i + 1] <<" has I_Q_2 density "<< y[i + 1 + num_parasites] << endl;
                //     cout << " and I_F_1 density " << y[i + 1 +2*num_parasites] << " and I_F_2 density " << y[i + 1 + 3*num_parasites] << endl;
                        
                // }
                // cout << "With " << y[3*num_parasites + 1] << " recovered hosts" << endl;
                if(y[0] <= TINY){
                    y[0] = 0.0;
                }
                for (i=0;i<num_parasites;i++){
                    float temp_sum = 0.0;
                    temp_sum += y[i+1];
                    temp_sum += y[1 + num_parasites + i];
                    temp_sum += y[1 + 2*num_parasites + i];
                    if (temp_sum <= TINY){
                        y[i+1] = 0.0;
                        y[1 + num_parasites + i] = 0.0;
                        y[1 + 2*num_parasites + i] = 0.0;
                    }
                }
                // for (i=0; i<(2 + 3*num_parasites ); i++){
                //     if(y[i] <= TINY){
                //         y[i] = 0.0;
                //     }
                // }
                t[0]+=h[0];
                h[0] = hnext[0];
                
                for(i=0;i<(2 + 3*num_parasites );i++){
                    y_max[i] = fastmax(y[i], y_max[i]);
                    y_min[i] = fastmin(y[i], y_min[i]);
                }
                step_count += 1;
                
                if ((((t[0] - tcheck[0]) >= 10.0) && (step_count - check_count) >= 200) && (y[0] != 0.0)){
                    discrep_check = 0.0;
                    check_count = step_count;
                    for(i=0;i<(2 + 3*num_parasites );i++){
                        if (y[i] >= TINY){
                            if (discrep_check < fastmax(abs(y_max[i] - y_max_next[i]),abs(y_min[i] - y_min_next[i]))){
                                discrep_check = fastmax(abs(y_max[i] - y_max_next[i]),abs(y_min[i] - y_min_next[i]));
                            }
                        }
                    }

                    if (discrep_check >= 1e-3){
                        y_check_sum = 0.0;

                        tcheck[0] = t[0];
                        for(i=0;i<(2 + 3*num_parasites );i++){
                            y_check[i] = y[i];
                            y_max[i] = y[i];
                            y_min[i] = y[i];
                            
                            y_max_next[i] = y[i];
                            y_min_next[i] = y[i];
                        }
                    }
                    else{
                        t[0] = tmax + 1.0;
                    }
                }
            }


            // cout << "After dynamics " << endl;
            // cout << y[0] << endl;
            // for(i = 0;i<num_parasites;i++){
            //     cout << "Parasite " << i <<" has I_Q_1 density "<< y[i + 1] <<" has I_Q_2 density "<< y[i + 1 + num_parasites] << endl;
            //     cout << " and I_F_1 density " << y[i + 1 +2*num_parasites] << " and I_F_2 density " << y[i + 1 + 3*num_parasites] << endl;
                        
            // }

            // cout << "With " << y[3*num_parasites + 1] << " recovered hosts" << endl;
            if (y[0] + y[3*num_parasites + 1] <= TINY){
                std::cout << "Healthy hosts extinct, breaking" << endl;
                break;
            }
            num_parasites_cleaned = 0;
            int para_counter = 0;
            int para_to_dos[trait_space_length];
            for(i=0;i<num_parasites;i++){
                float temp_sum = 0.0;
                temp_sum += y[i+1];
                temp_sum += y[1 + num_parasites + i];
                temp_sum += y[1 + 2*num_parasites + i];

                if (temp_sum > TINY){
                    num_parasites_cleaned += 1;
                    rho_inds_cleaned[para_counter] = rho_inds[i];
                    Q_temp[para_counter] = y[i + 1];
                    A_temp[para_counter] = y[i + 1 + num_parasites];
                    U_temp[para_counter] = y[i + 1 + 2*num_parasites];
                    para_to_dos[para_counter] = i;
                    para_counter++;
                }
            }

            if (num_parasites_cleaned == 0){
                cout << "All parasites are dead\n";
                break;
            }

            float rand_stored1, rand_stored2;

            
            for (i = 0; i<num_parasites_cleaned; i++){
                para_density[i] = Q_temp[i] + U_temp[i];
                para_density[i] += A_temp[i];
                if (i == 0){
                    para_cum_density[i] = para_density[i];
                }
                else{
                    para_cum_density[i] = para_cum_density[i-1] + para_density[i];
                }
            }

            for (i = 0; i< num_parasites_cleaned; i++){
                para_props[i] = para_cum_density[i]/para_cum_density[num_parasites_cleaned-1];
            }

            rand_stored1 = (float(rand())/float((RAND_MAX)));
            rand_stored2 = (float(rand())/float((RAND_MAX)));

            bool flag = true;
            int increment_ind;

            num_parasites_2 = num_parasites_cleaned;
            

            // cout << "What are the y values\n";
            // for(i = 0; i<(2 + 3*num_parasites ); i++){
            //     cout << "y_" << i << " = " << y[i] << endl;
            // }
            // cout << "Para check is\n";
            // cout << para_check << endl;

            for(i = 0; i<num_parasites_cleaned; i++){
                if (para_props[i] >= rand_stored1){
                    ind = i;
                    break;
                }
            }
            if ((((rand_stored2) >= 0.5) && (rho_inds_cleaned[ind] != (trait_space_length - 1))) || (rho_inds_cleaned[ind] == 0)){
                // cout << "Parasite going up" << endl;
                int ind_to_check = rho_inds_cleaned[ind] + 1;
                for(i = 0; i < num_parasites_cleaned; i++){
                    if (rho_inds_cleaned[i] == ind_to_check){
                        flag = false;
                        increment_ind = i;
                        break;
                    }
                }
            }
            else{
                // cout << "Parasite going down" << endl;
                int ind_to_check = rho_inds_cleaned[ind] - 1;
                for(i = 0; i < num_parasites_cleaned; i++){
                    if (rho_inds_cleaned[i] == ind_to_check){
                        flag = false;
                        increment_ind = i;
                        break;
                    }
                }
            }
            if(flag){
                num_parasites_2++;
            }

            y_temp[0] = y[0];
            for (i=1; i<(2 + 3*num_parasites_2); i++){
                y_temp[i] = 0.0;
            }
            y_temp[1 + 3*num_parasites_2] = y[1 + 3*num_parasites];
            // cout << "Filling the temp vector" << endl;
            // cout << y_temp[0] << endl;
            for (i=0; i<num_parasites_cleaned; i++){
                y_temp[i + 1] = Q_temp[i];
                y_temp[i + num_parasites_2 + 1] = A_temp[i];
                y_temp[i + 2*num_parasites_2 + 1] = U_temp[i];
                rho_inds_2[i] = rho_inds_cleaned[i];
                
            }
            num_parasites = num_parasites_2;
            
            for(i = 0; i<(2 + 3*num_parasites ); i++){
                y[i] = 0.0;
            }

            y[0] = y_temp[0];            
            y[1 + 3*num_parasites] = y_temp[1 + 3*num_parasites_2];

            for(i=0;i<num_parasites;i++){
                rho_inds[i] = rho_inds_2[i];
                y[i+1] = y_temp[i+1];
                y[i + num_parasites + 1] = y_temp[i + num_parasites + 1];
                y[i + 2*num_parasites + 1] = y_temp[i + 2*num_parasites + 1];
            }


            if ((((rand_stored2) >= 0.5) && (rho_inds_cleaned[ind] != (trait_space_length - 1))) || (rho_inds_cleaned[ind] == 0)){
                if(flag){
                    rho_inds[num_parasites-1] = rho_inds_2[ind] + 1;
                    y[num_parasites] = (y_temp[ind+1]/100);
                    y[2*num_parasites] = (y_temp[num_parasites+ind+1]/100);
                    y[3*num_parasites] = (y_temp[2*num_parasites+ind+1]/100);
                    
                }
                else{
                    y[1+increment_ind] = y[1+increment_ind] + (y_temp[ind+1]/100);
                    y[1+increment_ind + num_parasites] = y[1+increment_ind + num_parasites] + (y_temp[ind+1 + num_parasites]/100);
                    y[1+increment_ind + 2*num_parasites] = y[1+increment_ind + 2*num_parasites] + (y_temp[ind+1 + 2*num_parasites]/100);

                }
            }
            else{
                if(flag){
                    rho_inds[num_parasites-1] = rho_inds_2[ind] - 1;
                    

                    y[num_parasites] = (y_temp[ind+1]/100);
                    y[2*num_parasites] = (y_temp[num_parasites+ind+1]/100);
                    y[3*num_parasites] = (y_temp[2*num_parasites+ind+1]/100);
                    
                }
                else{
                    y[1+increment_ind] = y[1+increment_ind] + (y_temp[ind+1]/100);
                    y[1+increment_ind + num_parasites] = y[1+increment_ind + num_parasites] + (y_temp[ind+1 + num_parasites]/100);
                    y[1+increment_ind + 2*num_parasites] = y[1+increment_ind + 2*num_parasites] + (y_temp[ind+1 + 2*num_parasites]/100);

                }
            }
            
        
           
            for(i=0;i<num_parasites;i++){
                para_output_rho_inds[para_output_counter] = rho_inds[i];
                para_output_Q_densities[para_output_counter] = y[i + 1];
                para_output_A_densities[para_output_counter] = y[i + 1 + num_parasites];
                para_output_U_densities[para_output_counter] = y[i + 1 + 2*num_parasites];
                para_output_evosteps[para_output_counter] = evo_counter + 1;
                para_output_rho_vals[para_output_counter] = rho[rho_inds[i]];
                para_output_R_densities[para_output_counter] = y[3*num_parasites + 1];
                // para_output_beta_vals[para_output_counter] = beta[rho_inds[i]];
                para_output_host_density[para_output_counter] = y[0];
                para_output_counter++;
            }
            
            // cout << evo_counter + 1 << endl;
        }
        
        string test;

        string seed_str;
        seed_str = std::to_string(seed);

        // tracker_file.open("../data/tracker_file.csv", std::ios_base::app);
        // // tracker_file << "File_key,b,q,d,rho,sigma,gamma,eta,lambda"
        // tracker_file << "/coevo/data_set" + seed_str + ".csv," << b << "," << q << "," << d << "," << rho << "," << sigma_max << "," << gamma;
        // tracker_file << "," << eta << "," << lambda << "," << c1 << "," << c2 << "," << beta_max << "," << (hyper > 0.0) << "," << alpha_init << "," << sigma_init  << "\n";
        // tracker_file.close();

        // logs.close();

        
        para_outputs.open("../data/parasites/data_set"+seed_str+".csv");
        para_outputs << "Trait_index,host_density,density_of_Q,density_of_A,density_of_U,density_of_recovered_hosts,rho_value,beta_value,Evolutionary_step\n";
        for(i=0;i<para_output_counter;i++){
            para_outputs << para_output_rho_inds[i] <<"," << para_output_host_density[i] << "," << para_output_Q_densities[i] << "," << para_output_A_densities[i] << "," << para_output_U_densities[i] << "," << para_output_R_densities[i] << "," << para_output_rho_vals[i] << ',' << beta[para_output_rho_inds[i]] << "," << para_output_evosteps[i] << "\n";
        }
        para_outputs.close();
        
        delete y_check;
        delete y;
        delete yScale;
        delete dydt;
        delete y_temp;

        return 0;
    }

    void Quick_solver::eco_dynamics(float* S, float* Q, float* A, float* U, float*R, float beta_value, float rho_value, float b, float q, float d, float alpha, float gamma, float c1, float c2, float delta_Q, float eta, float sigma, int seed, float S_density, float Q_density, float A_density, float U_density){

        float y[5], yScale[5], dydt[5], discrep_check;
        float t[1], tcheck[1], tmax = 4000, h[1], hnext[1];
        float y_max[5], y_min[5], y_max_next[5], y_min_next[5], y_check_sum;
        int rho_inds[1];
        int num_parasites = 1;
        
        float beta[1], rho[1];
        int i, step_count = 0;
        rho_inds[0] = 0;
        

        y[0] = S_density/(b/q);
        y[1] = Q_density/(b/q);
        y[2] = A_density/(b/q);
        y[3] = U_density/(b/q);
        y[4] = 0.0;
        bool flag = true;

            

        rho[0] = rho_value;
        
        beta[0] = beta_value;
        

        rho_inds[0] = 0;
        
        
        t[0] = 0.0;
        h[0] = 1e-2;
        tcheck[0] = 10;
        for(i=0;i<(3*num_parasites + 2);i++){

            y_max[i] = y[i];
            y_min[i] = y[i];
            y_max_next[i] = y[i];
            y_min_next[i] = y[i];
            y_check_sum += y[i];
            // y_max[i] = y[i];
            // y_min[i] = y[i];
        }
        hnext[0] = 1e-2;

        while (t[0] <= tmax){   
            int check_count = 0;

            h[0] = hnext[0];
            
            /* This is where the equations are first solved */
            dynamics(dydt, y, num_parasites, rho_inds, b, q, d, beta, alpha, gamma, rho, delta_Q, eta, sigma);
        
            /* Adjust the step size to maintain accuracy */
            for (i=0; i<(3*num_parasites + 2); i++){
                yScale[i]=fabs(y[i])+fabs(dydt[i]*(*h))+TINY;
            }
            
            rkqs(y, dydt, h,hnext,yScale, num_parasites, rho_inds, b, q, d, beta, alpha, gamma, rho, delta_Q, eta, sigma, t);
            
            // cout << "During dynamics " << endl;
            // cout << "Host has density " <<  y[0] << " and" <<  endl;
            // for(i = 0;i<num_parasites;i++){
            //     cout << "Parasite " << i <<" has I_Q_1 density "<< y[i + 1] <<" has I_Q_2 density "<< y[i + 1 + num_parasites] << endl;
            //     cout << " and I_F_1 density " << y[i + 1 +2*num_parasites] << " and I_F_2 density " << y[i + 1 + 3*num_parasites] << endl;
                    
            // }

            // if (y[1] == 0.0 && y[2] == 0.0 && flag){
            //     flag = false;
            //     y[1] = TINY*10.0;
            // }

            t[0]+=h[0];

            for(i=0;i<(3*num_parasites + 2);i++){
                
                y_max[i] = fastmax(y[i], y_max[i]);
                
                y_min[i] = fastmin(y[i], y_min[i]);
                
            }
            // logs << "\n";
            step_count += 1;  
            
            if ((((t[0] - tcheck[0]) >= 10.0) || (step_count - check_count) >= 200) && (y[0] != 0.0)){
                discrep_check = 0.0;
                check_count = step_count;
                bool check_flag = false;
                
                for(i=0;i<(3*num_parasites + 2);i++){
                    if (discrep_check < fastmax(abs(y_max[i] - y_max_next[i]),abs(y_min[i] - y_min_next[i]))){
                        discrep_check = fastmax(abs(y_max[i] - y_max_next[i]),abs(y_min[i] - y_min_next[i]));
                        if (y_max[i] != 0.0 && y[i] == 0.0){
                            check_flag = true;
                        }
                    }
                }

                if (discrep_check >= TINY || check_flag){
                    y_check_sum = 0.0;

                    tcheck[0] = t[0];
                    for(i=0;i<(3*num_parasites + 2);i++){
                        y_max[i] = y[i];
                        y_min[i] = y[i];
                        
                        y_max_next[i] = y[i];
                        y_min_next[i] = y[i];
                    }
                }
                else{
                    t[0] = tmax + 1.0;;
                }
            }
        }

        S[0] = y[0];
        Q[0] = y[1];
        A[0] = y[2];
        U[0] = y[3];
        R[0] = y[4];
    }
}