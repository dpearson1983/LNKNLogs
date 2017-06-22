#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <random>
#include <cmath>
#include <fftw3.h>
#include <gsl/gsl_spline.h>
#include "../include/tpods.h"
#include "../include/constants.h"

// Return a vector with the FFT frequencies based on the grid properties
std::vector<double> fft_freq(int N, double L) {
    std::vector<double> k; // Declare this way with reserve for potential speed up
    k.reserve(N);
    double dk = (2.0*pi)/L; // Fundamental frequency
    for (int i = 0; i <= N/2; ++i)
        k.push_back(i*dk);
    for (int i = N/2 + 1; i < N; ++i)
        k.push_back((i - N)*dk);
    return k;
}

void fill_initial_Pk_grid(std::string in_pk_file, vec3<int> N, vec3<double> L, std::vector<double> &kx, 
                          std::vector<double> &ky, std::vector<double> &kz, std::vector<fftw_complex> dk, 
                          double b, double f) {
    std::ifstream fin;
    
    double V = L.x*L.y*L.z;
    
    std::vector<double> kin;
    std::vector<double> pin;
    
    if (std::ifstream(in_pk_file)) {
        fin.open(in_pk_file.c_str(), std::ios::in);
        while (!fin.eof()) {
            double k, p;
            fin >> k >> p;
            if (!fin.eof()) {
                kin.push_back(k);
                pin.push_back(p);
            }
        }
        fin.close();
    } else {
        std::stringstream message;
        message << "Could not open in_pk_file" << std::endl;
        throw std::runtime_error(message.str());
    }
    
    gsl_spline *Pk = gsl_spline_alloc(gsl_interp_cspline, pin.size());
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    
    gsl_spline_init(Pk, kin.data(), pin.data(), pin.size());
    
    for (int i = 0; i < N.x; ++i) {
        for (int j = 0; j < N.y; ++j) {
            for (int k = 0; k <= N.z/2; ++k) {
                double k_mag = sqrt(kx[i]*kx[i] + ky[j]*ky[j] + kz[k]*kz[k]);
                int index = k + (N.z/2 + 1)*(j + N.y*i);
                
                if (k_mag != 0) {
                    double mu = kz[k]/k_mag;
                    dk[index][0] = (b + mu*mu*f)*(b + mu*mu*f)*gsl_spline_eval(Pk, k_mag, acc)/V;
                    dk[index][1] = 0.0;
                } else {
                    dk[index][0] = 0.0;
                    dk[index][1] = 0.0;
                }
            }
        }
    }
    
    gsl_spline_free(Pk);
    gsl_interp_accel_free(acc);
}
