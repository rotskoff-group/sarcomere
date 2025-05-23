#ifndef MC_H
#define MC_H

#ifndef SARCOMERE_H
#include "sarcomere_old.h"
#endif

#ifndef UTILS_H
#include "utils.h"
#endif

#ifndef H5_UTILS_H
#include "h5_utils.h"
#endif

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <iostream>




class MC{
    public: 
        MC() = delete;
        MC(Sarcomere& model0, double& beta0, double& dt0, double& D0, int& update_myosin_every0,
            int& update_dt_every0, int& save_every0, bool& resume):
            model(model0), beta(beta0), dt(dt0), D(D0), update_myosin_every(update_myosin_every0), update_dt_every(update_dt_every0), save_every(save_every0){
            acc_rate=0;
            if (resume){
                model.load_state();
                printf("Resuming from file %s\n", model.filename.c_str());
                printf("Total energy: %f\n", model.get_energy());
                //print positions of myosins
                for (int i=0; i<model.myosin.n; i++){
                    printf("Myosin %d: %f %f\n", i, model.myosin.xs[i][0], model.myosin.xs[i][1]);
                    printf("Myosin endpoints: (%f %f), (%f %f)\n", model.myosin.left_endpts[i][0], model.myosin.left_endpts[i][1], model.myosin.right_endpts[i][0], model.myosin.right_endpts[i][1]);
                }
                //print positions of actins
                for (int i=0; i<model.actin.n; i++){
                    printf("Actin %d: %f %f\n", i, model.actin.xs[i][0], model.actin.xs[i][1]);
                }
                //print positions of alpha-actinins
                for (int i=0; i<model.alpha_actinin.n; i++){
                    printf("Alpha-actinin %d: %f %f\n", i, model.alpha_actinin.xs[i][0], model.alpha_actinin.xs[i][1]);
                }
                exit(0);
            }
            else{
                model.new_file();
            }
        }
        double dt, beta, acc_rate, D;
        int update_myosin_every, update_dt_every, save_every;
        Sarcomere& model;
        ~MC(){
        }

        void equilibrate(int myosin_nsteps, int alpha_actinin_nsteps, double dt0, gsl_rng * rng){
            double noise_mag = sqrt(2*D*dt0);
            double p_accept, delta_energy;
            for (int i=0; i<myosin_nsteps; i++){
                // gnerate random number from -1 to 1
                double randn[3];
                for (int i=0; i<3; i++)
                    randn[i]=gsl_rng_uniform(rng)*2-1;
                double acc_rand=gsl_rng_uniform(rng);
                // update myosin
                int myosin_index = gsl_rng_uniform_int(rng, model.myosin.n);
                double dx = randn[0]*noise_mag;
                double dy = randn[1]*noise_mag;
                double dtheta = randn[2]*M_PI*noise_mag;
                delta_energy = 0;
                model.equilibrate_myosin(myosin_index, dx, dy, dtheta, delta_energy);
                p_accept = exp(-beta*delta_energy);
                if (acc_rand>=p_accept){
                    model.equilibrate_myosin(myosin_index, -dx, -dy, -dtheta, delta_energy);
                }
            }
            for (int i=0; i<alpha_actinin_nsteps; i++){
                // gnerate random number from -1 to 1
                double randn[2];
                for (int i=0; i<2; i++)
                    randn[i]=gsl_rng_uniform(rng)*2-1;
                double acc_rand=gsl_rng_uniform(rng);
                // update alpha-actinin
                int alpha_actinin_index = gsl_rng_uniform_int(rng, model.alpha_actinin.n);
                double dx = randn[0]*noise_mag;
                double dy = randn[1]*noise_mag;
                delta_energy = 0;
                model.equilibrate_alpha_actinin(alpha_actinin_index, dx, dy, delta_energy);
                p_accept = exp(-beta*delta_energy);
                if (acc_rand>=p_accept){
                    model.equilibrate_alpha_actinin(alpha_actinin_index, -dx, -dy, delta_energy);
                }
            }
            model.total_energy = model.get_energy();
            printf("Total energy after equilibration: %f\n", model.total_energy);
        }

        void run_mc(int nsteps, gsl_rng * rng){
            double acc = 0;
            for (int i=0; i<nsteps; i++){
                bool update_myosin = (update_myosin_every>0 && i%update_myosin_every==0);
                acc+=sample_step(dt, D, rng, update_myosin);
                if (i>0 && update_dt_every>0 && i%update_dt_every==0){
                    acc_rate = acc/update_dt_every;
                    adjust_stepsize();
                    acc = 0;
                }
                if (i%save_every==0){
                    std::cout << "Step " << i << " Acceptance rate: " << acc_rate
                    << "Energy: "<<model.total_energy <<"Stepsize"<<dt<< std::endl;
                    model.save_state();
                    model.check_energy();
                }
            }
        }
        double sample_step(double& dt, double& D, gsl_rng * rng, bool& update_myosin){
            double acc = 0;
            double noise_mag = sqrt(2*D*dt);
            // gnerate random number from -1 to 1
            double randn[8];
            for (int i=0; i<8; i++)
                randn[i]=gsl_rng_uniform(rng)*2-1;
            double acc_rand[3];
            for (int i=0; i<3; i++)
                acc_rand[i]=gsl_rng_uniform(rng);
            // update actin
            int actin_index = gsl_rng_uniform_int(rng, model.actin.n);
            double dx = randn[0]*noise_mag;
            double dy = randn[1]*noise_mag;
            double dtheta = randn[2]*M_PI*noise_mag/10;
            //double dtheta = 0;
            double delta_energy = 0;
            std::vector <double> force(2);
            model.displace_actin(actin_index, dx, dy, dtheta, delta_energy, force);
            double work = force[0]*dx + force[1]*dy;
            double p_accept = exp(-beta*(delta_energy-work));
            if (acc_rand[0]<p_accept){
                acc+=1;
            }
            else{
                model.displace_actin(actin_index, -dx, -dy, -dtheta, delta_energy, force);
            }
            //printf("Actin displacement: %f %f %f\n", dx, dy, dtheta);
            if (update_myosin){
                // update myosin
                int myosin_index = gsl_rng_uniform_int(rng, model.myosin.n);
                dx = randn[3]*noise_mag;
                dy = randn[4]*noise_mag;
                dtheta = randn[5]*M_PI*noise_mag/10;
                delta_energy = 0;
                model.displace_myosin(myosin_index, dx, dy, dtheta, delta_energy,force);
                //printf("Myosin displacement: %f %f %f\n", dx, dy, dtheta);
                work = force[0]*dx + force[1]*dy;
                printf("Work: %f\n", work);
                p_accept = exp(-beta*(delta_energy-work));
                if (acc_rand[1]<p_accept){
                    acc+=1;
                }
                else{
                    model.displace_myosin(myosin_index, -dx, -dy, -dtheta, delta_energy, force);
                    //printf("Myosin displacement rejected\n");
                }
            }
            // update alpha_actinin
            int alpha_actinin_index = gsl_rng_uniform_int(rng, model.alpha_actinin.n);
            dx = randn[6]*noise_mag;
            dy = randn[7]*noise_mag;
            delta_energy = 0;
            model.displace_alpha_actinin(alpha_actinin_index, dx, dy, delta_energy);
            //printf("Alpha-actinin displacement: %f %f\n", dx, dy);
            p_accept = exp(-beta*delta_energy);
            if (acc_rand[2]<p_accept){
                acc+=1;
            }
            else{
                model.displace_alpha_actinin(alpha_actinin_index, -dx, -dy, delta_energy);
            }
            return acc/3;
        }
        void adjust_stepsize(){
            if (acc_rate<0.2){
                dt*=0.1;
            }
            else if (acc_rate<0.33){
                dt*=acc_rate/0.33;
            }
            else if (acc_rate>0.55){
                dt*=acc_rate/0.55;
            }	
            dt = std::max(dt, 1e-4);	
            dt = std::min(dt, 0.1);	
        }    
};

#endif

