#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <random>
#include <limits>
#include <cmath>
#include <stdexcept>
#include <fftw3.h>
#include <gsl/gsl_spline.h>
#include "../include/tpods.h"
#include "../include/constants.h"

// Return a vector with the FFT frequencies based on the grid properties
std::vector<double> fft_freq(int N, double L) {
    std::vector<double> k; // Declare this way with reserve for potential speed up
    k.reserve(N);
    double dk = (2.0*pi)/L; // Fundamental frequency
    
    // For the first half, the frequencies are positive and increasing
    for (int i = 0; i <= N/2; ++i)
        k.push_back(i*dk);
    
    // For the second half, the frequencies are negative, starting  and -k_max and going towards zero
    for (int i = N/2 + 1; i < N; ++i)
        k.push_back((i - N)*dk);
    
    return k;
}

// This fills an initial grid with the input power spectrum with plane-parallel anisotropies in the
// z direction.
void fill_initial_Pk_grid(std::string in_pk_file, vec3<int> N, vec3<double> L, std::vector<double> &kx, 
                          std::vector<double> &ky, std::vector<double> &kz, fftw_complex *dk, 
                          double b, double f) {
    // Setup input filestream
    std::ifstream fin;
    
    double V = L.x*L.y*L.z; // Volume of the grid
    
    // Vectors to store the input power spectrum data.
    std::vector<double> kin;
    std::vector<double> pin;
    
    // Check if the matter power specturm file exists, and if so, open it and read in the data. Otherwise,
    // throw an error and exit to prevent memory over-run problems.
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
    
    // Setup a cubic spline to distribute the input matter power to the grid
    gsl_spline *Pk = gsl_spline_alloc(gsl_interp_cspline, pin.size());
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    
    gsl_spline_init(Pk, kin.data(), pin.data(), pin.size());
    
    // Loop over the grid to distribute the input matter power with bias and anisotropy added in.
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

void normalize_dk_initial(fftw_complex *dk_i, int N, int N_tot) {
    for (int i = 0; i < N; ++i) {
        if (dk_i[i][0] > 0) {
            dk_i[i][0] /= N_tot;
            dk_i[i][1] = 0.0;
        } else {
            dk_i[i][0] = 0.0;
            dk_i[i][1] = 0.0;
        }
    }
}

std::vector<double> get_dk_initial(vec3<int> N, vec3<double> L, std::string wisdom_file, int nthreads,
                    std::vector<double> &kx, std::vector<double> &ky, std::vector<double> &kz, 
                    double b, double f, std::string in_pk_file) {
    int N_pad = N.x*N.y*2*(N.z/2 + 1);
    int N_rft = N.x*N.y*(N.z/2 + 1);
    int N_tot = N.x*N.y*N.z;
    fftw_init_threads();
    std::vector<double> dk_i(N_pad);
    
    fftw_import_wisdom_from_filename(wisdom_file.c_str());
    fftw_plan_with_nthreads(nthreads);
    fftw_plan dk2dr = fftw_plan_dft_c2r_3d(N.x, N.y, N.z, (fftw_complex *)dk_i.data(), dk_i.data(),
                                            FFTW_MEASURE);
    fftw_plan dr2dk = fftw_plan_dft_r2c_3d(N.x, N.y, N.z, dk_i.data(), (fftw_complex *)dk_i.data(),
                                            FFTW_MEASURE);
    fftw_export_wisdom_to_filename(wisdom_file.c_str());
    
    fill_initial_Pk_grid(in_pk_file, N, L, kx, ky, kz, (fftw_complex *)dk_i.data(), b, f);
    
    fftw_execute(dk2dr);
    
    for (int i = 0; i < N_pad; ++i) {
        if (i < N_tot) dk_i[i] = log(1.0 + dk_i[i]);
        else dk_i[i] = 0.0;
    }
    
    fftw_execute(dr2dk);
    
    normalize_dk_initial((fftw_complex *)dk_i.data(), N_rft, N_tot);
    
    fftw_destroy_plan(dk2dr);
    fftw_destroy_plan(dr2dk);
    
    return dk_i;
}

// Fills a cube with a Gaussian random realizations of the initial lognormal field.
void get_dk_realization(fftw_complex *dk, fftw_complex *dk_i, 
                        std::vector<double> &kx, std::vector<double> &ky, std::vector<double> &kz,
                        vec3<int> N, vec3<double> L) {
    // Setup the needed random number generator. This uses the C++ random standard library which
    // requires C++11 or later.
    std::random_device seeder; // Random number generator that uses system entropy, use for seeding
    std::mt19937_64 gen(seeder()); // Mersenne Twister pseudo-random number generator
    std::normal_distribution<double> dist(0.0, 1.0); // Normal distribution centered on zero, with
                                                     // sigma = 1
    
    // Loop over the grid to get the random realization.
    for (int i = 0; i < N.x; ++i) {
        int i2 = (2*N.x - i) % N.x;
        for (int j = 0; j < N.y; ++j) {
            int j2 = (2*N.y - j) % N.y;
            for (int k = 0; k <= N.z/2; ++k) {
                int index1 = k + (N.z/2 + 1)*(j + N.y*i);
                
                // The following ensures that after the inverse FFT the resulting field is real
                if ((i == 0 || i == N.x/2) && (j == 0 || j == N.y/2) && (k == 0 || k == N.z/2)) {
                    dk[index1][0] = dist(gen)*sqrt(dk_i[index1][0]);
                    dk[index1][1] = 0.0;
                } else if (k == 0 || k == N.z/2) {
                    int index2 = k + (N.z/2 + 1)*(j2 + N.y*i2); // Point where k_1 = -k_2
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

// This function can be used to create a scaled version of the random realization for a tracer with 
// a different bias.
void scale_dk_realization(fftw_complex *dk_1, fftw_complex *dk_2, 
                          std::vector<double> &ratio, int N) {
    for (int i = 0; i < N; ++i) {
        dk_2[i][0] = ratio[i]*dk_1[i][0];
        dk_2[i][1] = ratio[i]*dk_1[i][1];
    }
}

// Create the catalog by Poisson sampling the real-space density field
void get_galaxies_from_dr(std::vector<double> &dr, vec3<int> N, vec3<double> L, double b, double nbar,
                          std::string out_file) {
    // Output file stream
    std::ofstream fout;
    
    // Setup some random number generator stuff
    std::random_device seeder;
    std::mt19937_64 gen(seeder());
    std::uniform_real_distribution<double> dist(0.0, 1.0); // Used for placing galaxies within cell
    
    // Calculate the grid spacing and the average number of galaxies per cell
    vec3<double> delr = {L.x/double(N.x), L.y/double(N.y), L.z/double(N.z)};
    double n = nbar*delr.x*delr.y*delr.z;
    
    // Find the mean and variance of the real-space field
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
    
    // Loop over the grid and Poisson sample to create the catalog
    fout.open(out_file.c_str(), std::ios::app); // Open output file in append mode for multi-tracer case
    fout.precision(std::numeric_limits<double>::digits10); // Output full precision
    for (int i = 0; i < N.x; ++i) {
        double x_min = i*delr.x;
        for (int j = 0; j < N.y; ++j) {
            double y_min = j*delr.y;
            for (int k = 0; k < N.z; ++k) {
                double z_min = k*delr.z;
                int index = k + N.z*(j + N.x*i);
                
                double density = n*exp(dr[index] - variance/2.0); // Exponentiate the field
                std::poisson_distribution<int> p_dist(density); // Poisson distribution to sample from
                int num_gals = p_dist(gen); // Poisson sampling
                
                // Output the galaxies in that cell, distributed in a uniform random way within the cell.
                for (int gal = 0; gal < num_gals; ++gal) {
                    fout << x_min + dist(gen)*delr.x << " " << y_min + dist(gen)*delr.y << " ";
                    fout << z_min + dist(gen)*delr.y << " " << b << " " << nbar << "\n";
                }
            }
        }
    }
    fout.close();
}
