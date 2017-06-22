#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <random>
#include <limits>
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
                          std::vector<double> &ky, std::vector<double> &kz, std::vector<fftw_complex> &dk, 
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

void get_dk_realization(std::vector<fftw_complex> &dk, std::vector<fftw_complex> &dk_i, 
                        std::vector<double> &kx, std::vector<double> &ky, std::vector<double> &kz,
                        vec3<int> N, vec3<double> L) {
    std::random_device seeder;
    std::mt19937_64 gen(seeder());
    std::normal_distribution<double> dist(0.0, 1.0);
    
    for (int i = 0; i < N.x; ++i) {
        int i2 = (2*N.x - i) % N.x;
        for (int j = 0; j < N.y; ++j) {
            int j2 = (2*N.y - j) % N.y;
            for (int k = 0; k <= N.z/2; ++k) {
                int index1 = k + (N.z/2 + 1)*(j + N.y*i);
                if ((i == 0 || i == N.x/2) && (j == 0 || j == N.y/2) && (k == 0 || k == N.z/2)) {
                    dk[index1][0] = dist(gen)*sqrt(dk_i[index1][0]);
                    dk[index1][1] = 0.0;
                } else if (k == 0 || k == N.z/2) {
                    int index2 = k + (N.z/2 + 1)*(j2 + N.y*i2);
                    dk[index1][0] = dist(gen)*sqrt(dk_i[index1][0]/2);
                    dk[index1][1] = dist(gen)*sqrt(dk_i[index1][0]/2);
                    
                    dk[index2][0] = dk[index1][0];
                    dk[index2][1] = -dk[index1][1];
                } else {
                    dk[index1][0] = dist(gen)*sqrt(dk_i[index1][0]/2);
                    dk[index1][1] = dist(gen)*sqrt(dk_i[index1][0]/2);
                }
            }
        }
    }
}

void scale_dk_realization(std::vector<fftw_complex> &dk_1, std::vector<fftw_complex> &dk_2, 
                          std::vector<double> &ratio) {
    int N = dk_1.size();
    for (int i = 0; i < N; ++i) {
        dk_2[i][0] = ratio[i]*dk_1[i][0];
        dk_2[i][1] = ratio[i]*dk_1[i][1];
    }
}

void get_galaxies_from_dr(std::vector<double> &dr, vec3<int> N, vec3<double> L, double b, double nbar,
                          std::string out_file) {
    std::ofstream fout;
    
    std::random_device seeder;
    std::mt19937_64 gen(seeder());
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    
    vec3<double> delr = {L.x/double(N.x), L.y/double(N.y), L.z/double(N.z)};
    double n = nbar*delr.x*delr.y*delr.z;
    
    int N_tot = dr.size();
    double mean = 0.0;
    for (int i = 0; i < N_tot; ++i)
        mean += dr[i];
    mean /= double(N_tot);
    
    double variance = 0.0;
    for (int i = 0; i < N_tot; ++i) {
        dr[i] -= mean;
        variance += (dr[i]*dr[i]);
    }
    variance /= double(N_tot - 1.0);
    
    fout.open(out_file.c_str(), std::ios::app);
    fout.precision(std::numeric_limits<double>::digits10);
    for (int i = 0; i < N.x; ++i) {
        double x_min = i*delr.x;
        for (int j = 0; j < N.y; ++j) {
            double y_min = j*delr.y;
            for (int k = 0; k < N.z; ++k) {
                double z_min = k*delr.z;
                int index = k + N.z*(j + N.x*i);
                
                double density = n*exp(dr[index] - variance/2.0);
                std::poisson_distribution<int> p_dist(density);
                int num_gals = p_dist(gen);
                
                for (int gal = 0; gal < num_gals; ++gal) {
                    fout << x_min + dist(gen)*delr.x << " " << y_min + dist(gen)*delr.y << " ";
                    fout << z_min + dist(gen)*delr.y << " " << b << " " << nbar << "\n";
                }
            }
        }
    }
    fout.close();
}
