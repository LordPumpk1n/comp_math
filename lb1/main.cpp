#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

const double a = 1.0;
const double b = 3.0;
const double c_val = 1.0;
const double d = 5.0;
const double r = 0.006;
const double s = 4.0;
const double x0 = -1.56;
const double I = 3.2;

double arr[39] = {0.0};
double temp[3] = {0.0};
double a_dopri5[] = {1.0/5.0,
    3.0/40.0, 9.0/40.0, 
    44.0/45.0,-56.0/15.0, 32.0/9.0, 
    19372.0/6561.0,-25360.0/2187.0, 64448.0/6561.0,-212.0/729.0, 
    9017.0/3168.0, -355.0/33.0, 46732.0/5247.0, 49.0/176.0, -5103.0/18656.0, 
    35.0/384.0, 0, 500.0/1113.0, 125.0/192.0, -2187.0/6784.0, 11.0/84.0,
    35.0/284.0, 0, 500.0/1113.0, 125.0/192.0, -2187.0/ 6784.0, 11.0/84.0, 0};

double a_dopri8[] = {0.056,
    0.021, 0.062,
    0.031, 0, 0.094,
    0.312, 0, -1.172, 1.172,
    0.037, 0, 0, 0.188, 0.150,
    0.048, 0, 0, 0.112, -0.026, 0.013,
    0.017, 0, 0, 0.388, 0.036, 0.197, -0.173,
    0.069, 0, 0, -0.634, -0.161, 0.139, 0.941, 0.212,
    0.184, 0, 0, -2.469, -0.291, -0.026, 2.848, 0.281, 0.124,
    -1.215, 0, 0, 16.673, 0.916, -6.057, -16.004, 14.849, -13.372, 5.134,
    0.259, 0, 0, -4.774, -0.435, -3.049, 5.578, 6.156, -5.062, 2.194, 0.135,
    0.822, 0, 0, -11.659, -0.758, 0.714, 12.076, -2.128, 1.990, -0.234, 0.176, 0, 
    0.042, 0, 0, 0, 0, -0.055, 0.239, 0.704, -0.760, 0.661, 0.158, -0.238, 0.25};

double c_dopri5[] = {0.2, 0.3, 0.8, 8.0/9.0, 1.0, 1.0};
double c_dopri8[] = {0.056, 0.083, 0.125, 0.312, 0.375, 0.148, 0.465, 0.565, 0.650, 0.925, 1.0, 1.0};

void hindmarshRose(double t, double* vars, double *dx) {
    double x = vars[0];
    double y = vars[1];
    double z = vars[2];
    
    dx[0] = y - a * x * x * x + b * x * x - z + I; // dx/dt
    dx[1] = c_val - d * x * x - y;                  // dy/dt
    dx[2] = r * (s * (x - x0) - z);                     // dz/dt
}

// RK4
void rk4_step(double t, double h, double* vars) {
    hindmarshRose(t, vars, &arr[0]);
    int k_index = 0;
    
    for (int i = 0; i < 3; i++) {
        temp[i] = vars[i] + 0.5 * h * arr[k_index];
        k_index++;
    }
    hindmarshRose(t + 0.5 * h, temp, &arr[k_index]); 
    
    for (int i = 0; i < 3; i++) {
        temp[i] = vars[i] + 0.5 * h * arr[k_index];
        k_index++;
    }
    hindmarshRose(t + 0.5 * h, temp, &arr[k_index]);
    
    for (int i = 0; i < 3; i++) {
        temp[i] = vars[i] + h * arr[k_index];
        k_index++;
    }
    hindmarshRose(t + h, temp, &arr[k_index]);
    
    for (int i = 0; i < 3; i++) {
        vars[i] += h * (arr[i] + 2 * arr[3 + i] + 2 * arr[6 + i] + arr[9 + i]) / 6.0;
    }
}

// DOPRI5
void dopri5_step(double t, double h, double *vars) {
    //vector<double> c{0.2, 0.3, 0.8, 8.0/9.0, 1.0, 1.0};
    
    hindmarshRose(t, vars, &arr[0]);
    
    int k_index = 0;
    int a_index = 0;

    for (int i = 0; i < 6; i++) { 
        temp[0] = vars[0];
        temp[1] = vars[1];
        temp[2] = vars[2];
        for (int k = 0; k <= i; k++) { 
            temp[0] += h * a_dopri5[a_index] * arr[k * 3];
            temp[1] += h * a_dopri5[a_index] * arr[k * 3 + 1];
            temp[2] += h * a_dopri5[a_index] * arr[k * 3 + 2];
            /*for (int j = 0; j < 3; j++) {
                temp[j] += h * a_dopri5[a_index] * arr[k * 3 + j];
            }*/
            a_index++;
        }
        hindmarshRose(t + c_dopri5[i] * h, temp, &arr[(i + 1) * 3]);
    }

    for (int i = 0; i < 7; i++) {
        vars[0] += h * a_dopri5[21 + i] * arr[i * 3];
        vars[1] += h * a_dopri5[21 + i] * arr[i * 3 + 1];
        vars[2] += h * a_dopri5[21 + i] * arr[i * 3 + 2];
        /*for (int j = 0; j < 3; j++) {
            vars[j] += h * a_dopri5[21 + i] * arr[i * 3 + j];
        }*/
    }

}

// DOPRI8
void dopri8_step(double t, double h, double* vars) {
    //vector<double> c {0.056, 0.083, 0.125, 0.312, 0.375, 0.148, 0.465, 0.565, 0.650, 0.925, 1.0, 1.0};
    hindmarshRose(t, vars, &arr[0]);
        
    double prod = 0;
    int count = 0;
    int k_index = 0;
    int a_index = 0;

    for (int i = 0; i < 12; i++) {
        temp[0] = vars[0];
        temp[1] = vars[1];
        temp[2] = vars[2];
        for (int k = 0; k <= i; k++) { // Исправлено: начинаем с k=0
            temp[0] += h * a_dopri8[a_index] * arr[k * 3];
            temp[1] += h * a_dopri8[a_index] * arr[k * 3 + 1];
            temp[2] += h * a_dopri8[a_index] * arr[k * 3 + 2];
            /*for (int j = 0; j < 3; j++) {
                temp[j] += h * a_dopri8[a_index] * arr[k * 3 + j];
            }*/
            a_index++;
        }
        hindmarshRose(t + c_dopri8[i] * h, temp, &arr[(i + 1) * 3]);
    }

    for (int i = 0; i < 13; i++) {
        vars[0] += h * a_dopri8[78 + i] * arr[i * 3];
        vars[1] += h * a_dopri8[78 + i] * arr[i * 3 + 1];
        vars[2] += h * a_dopri8[78 + i] * arr[i * 3 + 2];
        /*for (int j = 0; j < 3; j++) {
            vars[j] += h * a_dopri8[78 + i] * arr[i * 3 + j];
        }*/
    }

    
}

void func(void method(double, double, double*), const string& filename, double h, double tmax) {
    ofstream output(filename);
    output << "Time,x,y,z\n";
    
    double vars[] = {0.0, 0.0, 0.0};
    
    double t = 0.0;
    
    while (t <= tmax) {
        output << t << "," << vars[0] << "," << vars[1] << "," << vars[2] << "\n";
        method(t, h, vars);
        t += h;
    }
    
    output.close();
}

int main() {
    double h = 0.005;
    double tmax = 1000.0;
    /*
    ofstream output("hindmarsh_rose_dopri5.csv");
    output << "Time,x,y,z\n";
    
    double vars[] = {0.0, 0.0, 0.0};
    
    double t = 0.0;
    
    while (t <= tmax) {
        output << t << "," << vars[0] << "," << vars[1] << "," << vars[2] << "\n";
        dopri5_step(t, h, vars);
        t += h;
       
    }
    
    output.close();
    */

    func(rk4_step, "hindmarsh_rose_rk4.csv", h, tmax);
    func(dopri5_step, "hindmarsh_rose_dopri5.csv", h, tmax);
    func(dopri8_step, "hindmarsh_rose_dopri8.csv", h, tmax);
    
    cout << "Complete" << endl;
    return 0;
}