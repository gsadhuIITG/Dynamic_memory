#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define EPS 1e-8
#define NX 1001     // spatial points
#define NMU 51     // phenotypic points
#define L 10.0      // spatial domain size
#define Lmu 1.0    // phenotypic domain size
#define T_final 60.0  // time steps
#define dt 0.0001   // time step

#define num_var_para 2
double temp_var[num_var_para]; // temp_var hold the value for initialization parameters provided at start of simulation

// function to break a string and separate numeral character (initializers) based on the delimiter

void tokenize(const char *s, const char *del) {
    char str_copy[256];
    strncpy(str_copy, s, sizeof(str_copy));
    str_copy[sizeof(str_copy) - 1] = '\0';  // Ensure null termination

    char *token = strtok(str_copy, del);
    int i = 0;
    while (token != NULL && i < num_var_para) {  // assumes max 10 parameters
        temp_var[i++] = atof(token);
        token = strtok(NULL, del);
    }
}

/*void tokenize(std::string s, std::string del = " ")
{
    printf("Hello");
    int start = 0;
    int end = s.find(del);
    int i = 0;
    while (end != -1) {
        temp_var[i] = std::stod(s.substr(start, end - start));
        start = end + del.size();
        end = s.find(del, start);
        i++;
    }
    temp_var[i] = std::stod(s.substr(start, end - start));
}*/
///////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[]) {

    tokenize(argv[1], "_"); /// to separate parameters from input string 

    int i, j, k;
    double dx = L / (NX - 1);
    double dmu = Lmu / (NMU - 1);
    double t0=0.0;
    FILE *fp,*fp1,*fp2,*fp3,*fp4,*fp5,*fp6,*fp7,*fp8, *fp9, *fp12, *fp13, *fp15,*fp16,*fp17;




    // Parameters
    double beta;
    double D_c = 0.4, D_h = 0.01, D_mu = 0.0;
    double hapt_h =0.001;
    double lambda_n = 0.1, K = 1.0;
    double mu_nh = 0.5, gamma_n = 0.01, gamma_h = 0.01; //gamma_n = 0.001, gamma_h = 0.001;
    double delta_e = 15.0, lambda_e = 0.0; //ecm Parameters delta_e is degradation rate
    double Theta_n = 1.5, Theta_h = 0.5; // oxygen consumption rate
    double alpha = 0.8; // hypoxic memory terms
    double c_H = 0.25, c_N = 0.2;  // oxygen threshold values
    double c_inf = 1.0, mu_hn0 = 0.5;  // oxygen at boundary and phenotype
    double beta_1=temp_var[0], beta_2=temp_var[1];  //beta_1 in normoxic and beta_2 is hypoxic

    printf("Running analysis for beta1 %.2f, beta2 %.2f \n" , beta_1, beta_2);

    // Allocate arrays
    double n[NX], n_new[NX];
    double c[NX], c_new[NX];
    double e[NX], e_new[NX];
    double h[NX][NMU], h_new[NX][NMU];
    double H_all[NX];

    // double output_time_points[11] = {0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0};
    // Deciding number of times at which output will be saved
    double time_interval = 1.0;
    int num_time_outputs = T_final/time_interval+1;
    double output_time_points[num_time_outputs];
    for (i = 0; i <num_time_outputs; i++){
        output_time_points[i] = i*time_interval;
    }
    

    // Initial conditions
    for (i = 0; i < NX; i++) {
        //double x = i * dx;
        // e[i] = (x > 2) ? 0.5 : 0.0;
        c[i] = c_inf;
        for (j = 0; j < NMU; j++) h[i][j] = 0.0;
    }

    for (i=0;i<NX;i++){
        double x = i * dx;
        if(x<4.0){
            n[i] = 0.0;
            e[i]=0.5;
        }
        else if(x>4.0 && x<6.0){
            n[i] = 0.2;
            e[i]=0.0;
        }
        else{
            n[i] = 0.0;
            e[i]=0.5;
        }

    }

    char filename[100];

    /*sprintf(filename, "h_density_beta_%.2f_period_%.2f_at_time_%.2f.txt",beta,T,t0);
    fp = fopen(filename, "w");
    for (i = 0; i < NX; i++){
        for(j=0;j<NMU;j++){
            fprintf(fp, "%lf\t%lf\t%lf\n",i*dx,j*dmu, h[i][j]);
        }
    }
    fclose(fp);*/

    sprintf(filename, "output/velocity_beta1_%.2f_beta2_%.2f.txt",beta_1,beta_2);
    fp15=fopen(filename,"w");

    sprintf(filename, "output/tumor_volume_beta1_%.2f_beta2_%.2f.txt",beta_1,beta_2);
    fp12=fopen(filename,"w");

    sprintf(filename,"output/tumor_width_beta1_%.2f_beta2_%.2f.txt",beta_1,beta_2);
    fp13=fopen(filename,"w");

    double x_old_min=0.6;
    int time_indx = 0; // to track the time points to consider for the output
    k = 0; // tracking number of the iteration/time steps
    while(t0<T_final){


        // Compute total hypoxic density and f term
        double H_total[NX], f[NX][NMU];
        for (i = 0; i < NX; i++) {
            H_total[i] = 0.0;
            for (j = 0; j < NMU; j++) {
                H_total[i] += h[i][j];

            }
        }
        for (i = 0; i < NX; i++) {
            for (j = 0; j < NMU; j++) {
                double mu_val = j * dmu;
                if (-(alpha * (c_inf - c[i]) + mu_val - mu_hn0)<0.0){
                    beta=beta_2;
                }
                else{
                    beta=beta_1;
                }
                f[i][j] = -(alpha * (c_inf - c[i]) + mu_val - mu_hn0) / beta;
                if((j == 1 || j == 0) && f[i][1]<0.0){
                    f[i][j] = 0.0;
                }

            }
        }
        //printf("%lf\n", f[500][1]);

        for(i=0;i<NX;i++){
            H_all[i]=H_total[i];
        }

        double total_mass = 0.0;
        for(i=0;i<NX;i++){
            total_mass += (H_total[i]+n[i])*dx;
        }

        double n_total_mass = 0.0;
        for(i=0;i<NX;i++){
            n_total_mass += (n[i])*dx;
        }

        double h_total_mass = 0.0;
        for(i=0;i<NX;i++){
            h_total_mass += (H_total[i])*dx;
        }


        // Update normoxic cells
        for (i = 0; i < NX; i++) {
            double Hc =0.5*(1.0+tanh(100.0*(c[i]-c_H))); // (c[i] > c_H) ? 1.0 : 0.0;
            double Lc =0.5*(1.0+tanh(100.0*(c_H-c[i])));// (c_H > c[i]) ? 1.0 : 0.0;
            double Nc = 0.5*(1.0+tanh(100.0*(c_N-c[i]))); //(c[i] < c_N) ? 1.0 : 0.0;
            double source = lambda_n * n[i] * (1.0 - (n[i] + H_total[i] + e[i]) / K) * Hc;
            double h2n = 0.0;
            for (j = 0; j < NMU; j++)
                // h2n += (j * dmu) * h[i][j] * dmu;
                h2n += (j * dmu) * h[i][j];

            h2n *= Hc;
            n_new[i] = n[i] + dt * (source + h2n - mu_nh * n[i] * Lc - gamma_n * n[i] * Nc);
            //n_new[i] = n[i] + dt * (- mu_nh * n[i] * Lc);


        }

        // Update hypoxic cells


        for (i = 1; i < NX-1; i++) {
            for (j = 1; j < NMU - 1; j++) {
                double Hc =0.5*(1.0+tanh(100.0*(c[i]-c_H)));             //(c[i] > c_H) ? 1.0 : 0.0;
                double Nc =0.5*(1.0+tanh(100.0*(c_N-c[i])));
                double Lc =0.5*(1.0+tanh(100.0*(c_H-c[i])));
                double mu_val = j * dmu;
                //printf("%lf\n", f[500][1]);

                // Spatial diffusion and heptotaxis
                double dh_dx2_hepto = (0.5/(dx*dx))*(D_h*(2.0-((n[i+1]+n[i])+(H_total[i+1]+H_total[i])+(e[i+1]+e[i]))/K)*(h[i+1][j]-h[i][j])-D_h*(2-(n[i-1]+n[i]+H_total[i-1]+H_total[i]+e[i-1]+e[i])/K)*(h[i][j]-h[i-1][j])-hapt_h*(h[i+1][j]+h[i][j])*(e[i+1]-e[i])+hapt_h*(h[i][j]+h[i-1][j])*(e[i]-e[i-1]));


                // Phenotypic advection and diffusion
                double dfdmu = 0.0;
                if(j == 1){
                    if(f[i][j+1] > 0.0 && f[i][j] > 0.0){
                        dfdmu =(f[i][j] * h[i][j]) / (dmu); // f[i][0] = 0;
                    }
                    else if(f[i][j+1] < 0.0 && f[i][j] == 0.0){
                        dfdmu =(f[i][j+1] * h[i][j+1] - f[i][j] * h[i][j]) / (dmu);
                    }
                }
                else {
                    if(f[i][j-1]<=0.0 && f[i][j]<=0.0 && f[i][j+1]<=0.0){
                        dfdmu=(f[i][j+1] * h[i][j+1] - f[i][j] * h[i][j]) / (dmu);
                    }
                    else if(f[i][j-1]>=0.0 && f[i][j]>=0.0 && f[i][j+1]>0.0){
                        dfdmu =(f[i][j] * h[i][j] - f[i][j-1] * h[i][j-1]) / (dmu);
                    }
                    else if(f[i][j-1]>=0.0 && f[i][j]>0.0 && f[i][j+1]<0.0){
                        dfdmu =(-f[i][j-1] * h[i][j-1]) / (dmu); // assuming f[i][j] = 0;
                    }
                    else if(f[i][j-1]>0.0 && f[i][j]<0.0 && f[i][j+1]<=0.0){
                        dfdmu =(f[i][j+1] * h[i][j+1]) / (dmu); // assuming f[i][j] = 0;
                    }
                }

                double d2h_dmu2 = (h[i][j+1] - 2*h[i][j] + h[i][j-1]) / (dmu*dmu);

                // Source terms
                double delta=0.0;

                if(mu_val==mu_hn0){
                    delta=1.0;
                    //printf("True\n");
                }

                //double delta = (fabs(mu_val - mu_hn0) <dmu) ? 1.0 : 0.0;
                //printf("mu_val=%lf\n", mu_val);

                double gain = mu_nh * n[i] * Lc * delta;

                //printf("gain=%lf\t%lf\n", mu_nh * n[i] * Lc,gain);


                h_new[i][j] = h[i][j] + dt * (dh_dx2_hepto+D_mu*d2h_dmu2-dfdmu)+ dt*(gain- mu_val*h[i][j]*Hc - gamma_h * h[i][j] * Nc);
                //h_new[i][j] = h[i][j] + dt *(mu_nh * n[i] * Lc * delta);

             //if(i==500&&j==1)
                 //printf("%lf\t%0.10lf\t%0.10lf\t%lf\t%0.10lf\t%0.10lf\t%0.10lf\n",t0,h_new[i][j], (f[i][j+1] * h[i][j+1] - f[i][j] * h[i][j]) / (dmu), f[i][j+1] * h[i][j+1], f[500][1] * h[500][1], h[i][j], f[500][1]);


            }
            //printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",t0,h_new[500][1], h_new[500][2],h_new[500][3],h_new[500][10],h_new[500][20]);
        }
        for (i = 1; i < NX - 1; i++) {
            double dc_dx2 = (c[i+1] - 2*c[i] + c[i-1]) / (dx*dx);
            c_new[i] = c[i] + dt * (D_c * dc_dx2 - Theta_n * c[i] * n[i] - Theta_h * c[i] * H_total[i]);

        }

        // Update ECM
        for (i = 0; i < NX; i++) {
            double loss = delta_e * e[i] * H_total[i];
            double growth = lambda_e * e[i] * (1.0 - (n[i] + H_total[i] + e[i]) / K);
            e_new[i] = e[i] + dt * (growth - loss);
        }

        // Neumann and Dirichlet BCs
        c_new[0] = c_inf;
        c_new[NX-1] = c_inf;


        for (j = 0; j < NMU; j++) {
            h_new[0][j] = h_new[1][j];
            h_new[NX-1][j] = h_new[NX-2][j];
        }
        /*for (i = 0; i < NX; i++) {
            h_new[i][0] = h_new[i][1];
            h_new[i][NMU-1] = h_new[i][NMU-2];
        }*/

        /*for (i = 0; i < NX; i++) {
                   h_new[i][0] = h_new[i][1]*(D_mu/(D_mu+dmu*f[i][0]));
                    h_new[i][NMU-1] = h_new[i][NMU-2]*(D_mu/(D_mu-dmu*f[i][0]));
                 }*/


        //printf("%0.10lf\t%0.10lf\n",mu_nh * n[i] * Lc)
        t0=t0+dt;
        //printf("t=%0.10lf\t%lf\t%lf\t%.10lf\t%.10lf\n",t0, n[1],total_mass,n_total_mass,h_total_mass);


        // Copy updated values
        for (i = 0; i < NX; i++) {
            n[i] = n_new[i];
            e[i] = e_new[i];
            c[i] = c_new[i];
            for (j = 0; j < NMU; j++)
               h[i][j] = h_new[i][j];
        }
        double  x_min;
        double velocity1;


        for(i=0;i<NX;++i){
            if(n[i]+H_all[i]>0.0000000001)
                //min=m[i];
                //  double x_right=i*dx;
                //if(double x_right>0.6){
                x_min=i*dx;
            //}

        }
        
        if(k%5 == 0){
            velocity1=(double)(x_min-x_old_min)/5; //1 is written for j%2000
            x_old_min=x_min;
        }

            if(t0 >= output_time_points[time_indx]){

                fprintf(fp15,"%0.10lf\t%0.10lf\t%0.10lf\t%0.10lf\n", t0,x_old_min, x_min,velocity1);


                // Compute total hypoxic density and f term

                for (i = 0; i < NX; i++) {
                    H_total[i] = 0.0;
                    for (j = 0; j < NMU; j++) {
                        H_total[i] += h[i][j];
                    }
                }

                total_mass = 0.0;
                for(i=0;i<NX;i++){
                    total_mass += (H_total[i]+n[i])*dx;
                }

                fprintf(fp12,"%0.10lf\t%0.10lf\n",t0,total_mass);

                double width;
                width=2.0+2.0*(x_min-6.0);

                fprintf(fp13,"%0.10lf\t%0.10lf\n",t0,width);

                sprintf(filename, "output/n_density_beta1_%.2f_beta2_%.2f_at_time_%.2f.txt",beta_1,beta_2,t0);
                fp1 = fopen(filename, "w");

                sprintf(filename, "output/oxygen_conc_beta1_%.2f_beta2_%.2f_at_time_%.2f.txt",beta_1,beta_2,t0);
                fp2 = fopen(filename, "w");

                sprintf(filename, "output/e_density_beta1_%.2f_beta2_%.2f_at_time_%.2f.txt",beta_1,beta_2,t0);
                fp3 = fopen(filename, "w");

                sprintf(filename, "output/h_spatial_density_beta1_%.2f_beta2_%.2f_at_time_%.2f.txt",beta_1,beta_2,t0);
                fp4 = fopen(filename, "w");

                for (i = 0; i < NX; i++){
                    fprintf(fp1, "%lf\t%lf\n",i*dx, n[i]);
                    fprintf(fp2, "%lf\t%lf\n",i*dx, c[i]);
                    fprintf(fp3, "%lf\t%lf\n",i*dx, e[i]);
                    fprintf(fp4, "%lf\t%lf\n",i*dx, H_all[i]);
                }
                fclose(fp1);
                fclose(fp2);
                fclose(fp3);
                fclose(fp4);

                sprintf(filename, "output/h_density_beta1_%.2f_beta2_%.2f_at_time_%.2f.txt",beta_1,beta_2,t0);
                fp5 = fopen(filename, "w");
                for (i = 0; i < NX; i++){
                    for(j=0;j<NMU;j++){
                        fprintf(fp5, "%lf\t%lf\t%lf\n",i*dx,j*dmu, h[i][j]);
                    }
                }
                fclose(fp5);

                time_indx++;
            }

        k++; // incrementing the count of number of time steps evaluated
    }

    fclose(fp15);
    fclose(fp12);
    fclose(fp13);
    return 0;
}
