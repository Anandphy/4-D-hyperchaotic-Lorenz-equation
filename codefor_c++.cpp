// Bifurcation diagram and time series analysis for 4-D Lorenz equation.
// R. Anand, anandphy0@gmail.com

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <array>

const int n = 4;
const int j_max = 10;
const int k_max = 100000;
const int j_store = 10;
const int n_run = 10000;
const int n_pas = 5000;
const double a = 0.87;
const double b = 0.5;
const double d = 1.5;
const double m = 0.5;
const double stp_size = 1.0 / j_max;
const double errlt = 1.0e-12;
const int neps = 1000;
const double c_min = 0.0;
const double c_max = 3.0;
const double c_stp = (c_max - c_min) / neps;

double c;

std::array<double, n> ode(const std::array<double, n>& z, double t) {
    std::array<double, n> res;
    res[0] = a * (z[1] - z[0]);
    res[1] = -b * z[1] + z[0] * z[2] + c * z[3];
    res[2] = d - z[0] * z[1];
    res[3] = -m * z[0];
    return res;
}

void runge_kutta(std::array<double, n>& x, double t, double h, std::array<double, n> (*f)(const std::array<double, n>&, double)) {
    std::array<double, n> k1, k2, k3, k4, tmp;
    k1 = f(x, t);
    for (int i = 0; i < n; ++i) {
        k1[i] *= h;
        tmp[i] = x[i] + k1[i] / 2.0;
    }
    k2 = f(tmp, t + h / 2.0);
    for (int i = 0; i < n; ++i) {
        k2[i] *= h;
        tmp[i] = x[i] + k2[i] / 2.0;
    }
    k3 = f(tmp, t + h / 2.0);
    for (int i = 0; i < n; ++i) {
        k3[i] *= h;
        tmp[i] = x[i] + k3[i];
    }
    k4 = f(tmp, t + h);
    for (int i = 0; i < n; ++i) {
        k4[i] *= h;
        x[i] += (k1[i] + 2.0 * (k2[i] + k3[i]) + k4[i]) / 6.0;
    }
}

void init_vals(std::array<double, n>& x, std::array<double, n>& y, double& dt, double& t) {
    dt = stp_size;
    x = {0.1, 0.1, 0.1, 0.1};
    for (int i = 0; i < n_pas; ++i) {
        for (int j = 0; j < j_max; ++j) {
            runge_kutta(x, t, dt, ode);
            t += dt;
        }
    }
    y = x;
}

int main() {
    std::ofstream time("ts.dat"), bif("bifurcation.txt");

    std::array<double, n> x, y;
    std::vector<double> peaks(k_max);
    double t, dt, t1;
    int k_count, ij;

    for (int i = 0; i <= neps; ++i) {
        c = c_min + i * c_stp;
        std::cout << c << std::endl;
        init_vals(x, y, dt, t);
        k_count = 0;
        for (int j = 1; j <= n_run; ++j) {
            for (ij = 1; ij <= j_max; ++ij) {
                runge_kutta(x, t, dt, ode);
                time << t << " " << x[0] << " " << x[1] << std::endl;
                //time2 << x[2] << " " << x[3] << std::endl;
                if (y[1] > 0.0 && x[1] < 0.0) {
                    if (std::abs(x[1]) < errlt) {
                        t += dt;
                        k_count++;
                        peaks[k_count] = x[0];
                        bif << x[0] << " " << c << std::endl;
                        dt = stp_size;
                        y = x;
                    } else {
                        x = y;
                        dt /= 2.0;
                    }
                } else {
                    y = x;
                    t += dt;
                }
            }
        }
    }

    time.close();
    //time2.close();
    bif.close();

    return 0;
}

