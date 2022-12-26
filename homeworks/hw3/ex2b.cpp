#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <random>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <format>
#include <random>
using namespace std;

vector <double> compute_F_i(double tau,int N, double A, double B, double C, double D, int i,vector <double> v_x0, vector <double> v_y0, vector <double> x, vector <double> y, vector <double> v_x, vector <double> v_y,  vector <double> o_x, vector <double> o_y);
double average_x_displacement(vector <double>x, vector <double>x0);
vector <double> differential_evolution(int NP, double CR, double F, double tau, int N, int max_time,
                                    double A, double B, double C, double D, vector <double> v_x0,
                                    vector <double> v_y0, vector <double> x, vector <double> y, vector <double> x0 ,vector <double> v_x, 
                                    vector <double> v_y,  vector <double> o_x, vector <double> o_y, 
                                    double delta_t, double high_y, double low_y, mt19937 gen);
double f(vector <double> agent, vector <double> x0,vector <double> x, vector <double> y, 
        vector <double> v_x0,vector <double> v_y0, vector <double> v_x, vector <double> v_y, 
        int N, int max_time,  double tau,  double A, double B, double C, double D, double delta_t, 
        double high_y, double low_y );
int main()
{
    //hyperparameters

    int N = 10;
    int NP = 8;
    double CR = 0.9;
    double F = 0.8;
    double tau = 0.2;
    double A = 20.;
    double B = 0.5;
    double C = 10.;
    double D = 0.6;
    double delta_t = 0.05;
    int max_time = 1000;
    double high_y = 3.;
    double low_y = -3;
    mt19937 gen(20);

   // initialise positions
    vector <double> x {-28.5, -27.0, -25.5, -24.0, -22.5, 22.5, 24.0, 25.5, 27.0, 28.5};
    vector <double> y {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};

    vector <double>  x0  {-28.5, -27.0, -25.5, -24.0, -22.5, 22.5, 24.0, 25.5, 27.0, 28.5};
    vector <double>  y0 {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};

    // initialise velocities
    vector <double> v_x0 {1., 1., 1., 1., 1., -1, -1, -1, -1, -1};
    vector <double> v_y0 {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};

    vector <double> v_x {1., 1., 1., 1., 1., -1, -1, -1, -1, -1};
    vector <double> v_y  {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};

    // create obstacles 
    vector <double> o_x (5);
    vector <double> o_y (5);
    
    uniform_real_distribution < double> distrib_x(-5 + D,5 - D); 
    uniform_real_distribution < double> distrib_y(-3 + D,3 - D); 

    for (int i = 0; i < 5 ; i++){
            o_x[i] = distrib_x(gen);
		}
    for (int i = 0; i < 5 ; i++){
            o_y[i] = distrib_y(gen);
		}

    vector <double> best_obstacles(10);
    best_obstacles = differential_evolution(NP, CR, F, tau, N, max_time, A, B, C, D, v_x0, v_y0, x, y,x0, v_x, v_y, o_x, o_y, delta_t, high_y, low_y, gen);
    
    cout << "\nbest obstacles : [";
    for(int i = 0;i < 5; i++){
        if(i ==4){
            cout << "["<<best_obstacles[i] << "," << best_obstacles[i + 5]<< "]";
        }
        else{
            cout << "["<<best_obstacles[i] << "," << best_obstacles[i + 5]<< "],\n";
        }
    }
    cout << "]\n\n";
    double best_average_displacement = f(best_obstacles, x0,x,y, v_x0,v_y0, v_x, v_y, N,  max_time, tau,A, B, C, D, delta_t, high_y, low_y );
    cout << "maximal average displacement is : " <<best_average_displacement<< "\n\n";

    cout << "obstacle for c++ : \n{";
    for(int i = 0;i < 10; i++){
        if(i <4){
            cout <<best_obstacles[i] << ",";
        }
        if(i == 4){
            cout <<best_obstacles[i] << "}\n{";
        }
        if (i >= 5 and i < 9){
            cout <<best_obstacles[i]<< ",";
        }
        if (i == 9){
            cout <<best_obstacles[i];
        }
    }
    cout << "}\n\n";


    return 0;
}

vector <double> compute_F_i(double tau, int N, double A, double B, double C, double D, int i,vector <double> v_x0, vector <double> v_y0, vector <double> x, vector <double> y, vector <double> v_x, vector <double> v_y,  vector <double> o_x, vector <double> o_y){
    
    vector <double> F_i (2);
    double rik = 0;
    double rij = 0;
    double sum_x = (1/tau)*(v_x0[i] - v_x[i]);
    double sum_y = (1/tau)*(v_y0[i] - v_y[i]);
    for (int j = 0; j< N; j++){
        if (i != j){
            rij = sqrt(pow(x[i] - x[j],2) + pow(y[i] - y[j], 2));
            sum_x += A*exp(-rij/B)*(x[i] - x[j])/rij;
            sum_y += A*exp(-rij/B)*(y[i] - y[j])/rij;
        }
    }

    for (int k = 0; k< 5; k++){
        rik = sqrt(pow(x[i] - o_x[k],2) + pow(y[i] - o_y[k], 2));
        sum_x += C*exp(-rik/D)*(x[i] - o_x[k])/rik;
        sum_y += C*exp(-rik/D)*(y[i] - o_y[k])/rik;
    }
    F_i[0] = sum_x;
    F_i[1] = sum_y;
    return F_i;
}

double average_x_displacement(vector <double>x, vector <double>x0){
    double sum = 0;
    for(int i = 0;i < 10; i++){
        sum += abs(x[i] - x0[i]);
    }
    return sum/10;
}

vector <double> differential_evolution(int NP, double CR, double F, double tau, int N, int max_time,
                                    double A, double B, double C, double D, vector <double> v_x0,
                                    vector <double> v_y0, vector <double> x, vector <double> y, vector <double> x0, vector <double> v_x, 
                                    vector <double> v_y,  vector <double> o_x, vector <double> o_y, 
                                    double delta_t, double high_y, double low_y, mt19937 gen){
    // intialise agents obstacles
 
    uniform_real_distribution < double> distrib_x(-5 + D,5 - D);
    uniform_real_distribution < double> distrib_y(-3 + D,3 - D); 
    uniform_real_distribution < double> distrib_r(0,1);
    uniform_int_distribution < int> distrib_agent(0,NP);
    uniform_int_distribution < int> distrib_R(0,10);
    
    double population[N][NP]; //cols, rows
    for (int coord = 0; coord < N; coord++){
        for (int pop = 0; pop < NP; pop++){
            if(coord < 5){
                population[coord][pop] = distrib_x(gen);
            }
            if(coord >=5){
                population[coord][pop] = distrib_y(gen);
            }
        }
    }

    // termination criterion
    int max_iter = 100;
    int t = 0;
    while(t < max_iter){
        for (int i = 0;i < NP; i++){ // for each agent from agent population
            vector <double> current_agent(10);

            for(int j = 0; j < 10; j++){
                current_agent[j] = population[j][i];
            }

            // make 3 agents a,  b, c
            vector <double> a(10);
            vector <double> b(10);
            vector <double> c(10);

            int a_idx = distrib_agent(gen);
            int b_idx = distrib_agent(gen);
            int c_idx = distrib_agent(gen);
            while(a_idx == b_idx or a_idx == c_idx or c_idx == b_idx or a_idx== i or b_idx == i or c_idx == i){
                b_idx = distrib_agent(gen);
                c_idx = distrib_agent(gen);
                a_idx = distrib_agent(gen);
            }

            for(int j = 0; j < 10; j++){
                a[j] = population[j][a_idx];
                b[j] = population[j][b_idx];
                c[j] = population[j][c_idx];
            }

            //get random index in [0,10]
            int R = distrib_R(gen);

            //new potential position for agent
            vector <double> y(10);

            for (int j =0; j < 10; j++){
                double r_j = distrib_r(gen);
                if(r_j < CR or j == R){
                    vector <double> b_c(10);
                    for (int k = 0; k < 10 ; k++){
                        b_c[k] = F*(b[k] - c[k]);
	                }
                    
                    y[j] = a[j]+b_c[j];
                    // cout << y[j]<< "\n";
                }
                if (r_j >= CR and j != R){
                    y[j] = current_agent[j];
                }
            }
            double distance_y = f(y, x0,x,y, v_x0,v_y0, v_x, v_y, N,  max_time, tau,A, B, C, D, delta_t, high_y, low_y );
            double distance_current_agent = f(current_agent, x0,x,y, v_x0,v_y0, v_x, v_y, N,  max_time, tau,A, B, C, D, delta_t, high_y, low_y );

            if ( distance_y> distance_current_agent){
                for(int j = 0; j < 10; j++){
                    if (j <5 and abs(y[j]) <= 5 - D ){
                        population[j][i] = y[j];
                    }
                    if (j >= 5 and abs(y[j])<= 3 - D ){
                        population[j][i] = y[j];
                    }
                    
                }
            }
        }
        t += 1;
    }

    double distance = -1;
    vector <double> best_agent(10);
    for (int i = 0; i < NP; i++ ){
        vector <double> agent(10);
        for(int j = 0; j < 10; j++){
            agent[j] = population[j][i];
        }
        double cur_dist = f(agent, x0, x,y, v_x0,v_y0, v_x, v_y, N,  max_time, tau,A, B, C, D, delta_t, high_y, low_y );
        if(  cur_dist > distance ){
            for(int j = 0; j < 10; j++){
                best_agent[j] = agent[j];
            }
        }
    }

    return best_agent;
}

double f(vector <double> agent, 
        vector <double> x0,
        vector <double> x, 
        vector <double> y, 
        vector <double> v_x0,
        vector <double> v_y0, 
        vector <double> v_x, 
        vector <double> v_y, 
        int N, 
        int max_time, 
        double tau, 
        double A, 
        double B, 
        double C, 
        double D, 
        double delta_t, 
        double high_y, 
        double low_y){

    // create obstacles 
    vector <double> o_x (5);
    vector <double> o_y (5);


    for (int i = 0; i < 5 ; i++){
            o_x[i] = agent[i];
		}
    for (int i = 0; i < 5 ; i++){
            o_y[i] = agent[i];
		}
 
    // initialise force

    vector <double> F_i (2);

    for (int time = 0; time < max_time ; time++){
        for(int i = 0;  i < N; i++){
            F_i   = compute_F_i(tau, N, A, B, C, D, i, v_x0, v_y0, x, y, v_x, v_y, o_x, o_y);
            double F_i_x = F_i[0];
            double F_i_y = F_i[1];
            x[i] += v_x[i]*delta_t;
            y[i] += v_y[i]*delta_t;
            v_x[i] += F_i_x*delta_t;
            v_y[i] += F_i_y*delta_t;

            if (y[i]>high_y){
                y[i] = high_y;
                }
            if (y[i]<low_y){
                y[i] = low_y;
            }
        }
    }
    double distance = 0;
    distance = average_x_displacement(x, x0);
    return distance;
}

