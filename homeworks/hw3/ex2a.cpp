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


int main()
{
    //hyperparameters

    int N = 10;
    double tau = 0.2;
    double A = 20.;
    double B = 0.5;
    double C = 10.;
    double D = 0.6;
    int O = 5;
    double delta_t = 0.05;
    int max_time = 1000;
    double high_y = 3.;
    double low_y = -3;
    mt19937 gen(20);

   // initialise positions
    vector <double> x {-28.5, -27.0, -25.5, -24.0, -22.5, 22.5, 24.0, 25.5, 27.0, 28.5};
    vector <double>  y {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};

    vector <double>  x0  {-28.5, -27.0, -25.5, -24.0, -22.5, 22.5, 24.0, 25.5, 27.0, 28.5};
    vector <double>  y0 {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};

    // initialise velocities
    vector <double> v_x0 {1., 1., 1., 1., 1., -1, -1, -1, -1, -1};
    vector <double> v_y0 {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};

    vector <double> v_x {1., 1., 1., 1., 1., -1, -1, -1, -1, -1};
    vector <double> v_y  {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};

    // create obstacles 

    // if you already know the positions of your obstacles
    // vector <double> o_x {-3.39921,0.0510252,-0.29937,-4.27294,2.00433};

    // vector <double> o_y {1.08641,2.04437,-0.361777,0.268635,-1.01766};


    // generate random obstacles
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

    cout << "Obtsacles : [";
    for(int i = 0;i < 5; i++){
        if(i == 4){
            cout << "["<<o_x[i] << "," << o_y[i]<< "]";
        }
        else{
            cout << "["<<o_x[i] << "," << o_y[i]<< "],\n";
        }
    }
    cout << "]\n\n";
    
    // initialise force

    vector <double> F_i (2);

    for (int time = 0; time < max_time ; time++){
        for(int i = 0;  i < N; i++){
            F_i   = compute_F_i(tau, N, A, B, C, D, i, v_x0, v_y0, x, y, v_x, v_y, o_x, o_y);
            // cout << "F_i_y is :" << F_i[1]<< "\n";
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
    cout << "Final Points:\n [";
    for(int i = 0;i < 10; i++){
        if(i == 9){
            cout << "["<<x[i] << "," << y[i]<< "]";
        }
        else{
            cout << "["<<x[i] << "," << y[i]<< "],\n";
        }
    }
    cout << "]\n";

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

