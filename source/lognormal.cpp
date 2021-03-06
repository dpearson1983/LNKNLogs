/* LNKNLogs v. 1.0
 * David W. Pearson & Lado Samushia
 * Copyright (C) 2016, 2017
 * 
 * Contains function definitions for LNKNLogs.
 * 
 * LICENSE: GPL v3
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 */

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
#include "../include/lognormal.h"
#include "../include/tpods.h"
#include "../include/constants.h"
#include "../include/file_io.h"
#include "../include/cosmology.h"

#ifndef PI
#define PI 3.1415926535897932384626433832795
#endif

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
// x direction.
void fill_initial_Pk_grid(std::string in_pk_file, vec3<int> N, vec3<double> L, std::vector<double> &kx, 
                          std::vector<double> &ky, std::vector<double> &kz, std::vector<fftw_complex> &dk, 
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
                    double mu = kx[i]/k_mag;
                    dk[index][0] = ((b + mu*mu*f)*(b + mu*mu*f)*gsl_spline_eval(Pk, k_mag, acc))/V;
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

void normalize_dk_initial(std::vector<fftw_complex> &dk_i, int N, int N_tot) {
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

std::vector<fftw_complex> get_dk_initial(vec3<int> N, vec3<double> L, std::string wisdom_file, int nthreads,
                    std::vector<double> &kx, std::vector<double> &ky, std::vector<double> &kz, 
                    double b, double f, std::string in_pk_file) {
    int N_pad = N.x*N.y*2*(N.z/2 + 1);
    int N_rft = N.x*N.y*(N.z/2 + 1);
    int N_tot = N.x*N.y*N.z;
    fftw_init_threads();
    std::vector<double> dr_i(N_tot);
    std::vector<fftw_complex> dk_i(N_rft);
    
    fftw_import_wisdom_from_filename(wisdom_file.c_str());
    fftw_plan_with_nthreads(nthreads);
    fftw_plan dk2dr = fftw_plan_dft_c2r_3d(N.x, N.y, N.z, dk_i.data(), dr_i.data(), FFTW_MEASURE);
    fftw_plan dr2dk = fftw_plan_dft_r2c_3d(N.x, N.y, N.z, dr_i.data(), dk_i.data(), FFTW_MEASURE);
    fftw_export_wisdom_to_filename(wisdom_file.c_str());
    
    fill_initial_Pk_grid(in_pk_file, N, L, kx, ky, kz, dk_i, b, f);
    
    fftw_execute(dk2dr);
    
    for (int i = 0; i < N_tot; ++i)
        dr_i[i] = log(1.0 + dr_i[i]);
    
    fftw_execute(dr2dk);
    
    normalize_dk_initial(dk_i, N_rft, N_tot);
    
    fftw_destroy_plan(dk2dr);
    fftw_destroy_plan(dr2dk);
    
    return dk_i;
}

// Fills a cube with a Gaussian random realizations of the initial lognormal field.
void get_dk_realization(std::vector<fftw_complex> &dk, std::vector<fftw_complex> &dk_i, 
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
                    dk[index1][0] = dist(gen)*sqrt(dk_i[index1][0]/2.0);
                    dk[index1][1] = dist(gen)*sqrt(dk_i[index1][0]/2.0);
                    
                    dk[index2][0] = dk[index1][0];
                    dk[index2][1] = -dk[index1][1];
                } else {
                    dk[index1][0] = dist(gen)*sqrt(dk_i[index1][0]/2.0);
                    dk[index1][1] = dist(gen)*sqrt(dk_i[index1][0]/2.0);
                }
            }
        }
    }
}

// This function can be used to create a scaled version of the random realization for a tracer with 
// a different bias.
void scale_dk_realization(std::vector<fftw_complex> &dk_1, std::vector<fftw_complex> &dk_2, 
                          std::vector<double> &ratio, int N) {
    for (int i = 0; i < N; ++i) {
        dk_2[i][0] = ratio[i]*dk_1[i][0];
        dk_2[i][1] = ratio[i]*dk_1[i][1];
    }
}

double arctangent(double x, double y) {
    if (x < 0 && y < 0) {
        return PI + atan(fabs(y)/fabs(x));
    } else if (x < 0) {
        return PI - atan(fabs(y)/fabs(x));
    } else if (y < 0) {
        return 2.0*PI - atan (fabs(y)/fabs(x));
    } else {
        return atan(fabs(y)/fabs(x));
    }
}

vec3<double> get_equatorial(vec3<double> cart, cosmology &cosmo) {
    double r = sqrt(cart.x*cart.x + cart.y*cart.y + cart.z*cart.z);
    vec3<double> equa;
    equa.y = atan2(cart.z, sqrt(cart.x*cart.x + cart.y*cart.y))*180.0/PI;
    equa.x = atan2(cart.y, cart.x)*180.0/PI;
    if (equa.x < 0) equa.x += 360.0;
    equa.z = cosmo.get_redshift_from_comoving_distance(r);
    return equa;
}

// Create the catalog by Poisson sampling the real-space density field
void get_galaxies_from_dr(std::vector<double> &dr, vec3<int> N, vec3<double> L, vec3<double> r_min,
                          double b, double nbar, cosmology &cosmo, std::string out_file, FileType type) {
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
    unsigned long N_tot = dr.size();
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
    unsigned long totalGals = 0;
    
    if (type == dr12) {
        std::vector<std::vector<double>> galaxies(8);
        std::vector<std::string> col_names;
        std::vector<std::string> col_forms;
        std::vector<std::string> col_units;
        get_DR12_column_info(col_names, col_forms, col_units);
        
        // Loop over the grid and Poisson sample to create the catalog
        //     fout.open(out_file.c_str(), std::ios::app); // Open output file in append mode for multi-tracer case
        //     fout.precision(std::numeric_limits<double>::digits10); // Output full precision
        for (int i = 0; i < N.x; ++i) {
            double x_min = i*delr.x;
            for (int j = 0; j < N.y; ++j) {
                double y_min = j*delr.y;
                for (int k = 0; k < N.z; ++k) {
                    double z_min = k*delr.z;
                    size_t index = k + N.z*(j + N.x*i);
                    
                    double density = n*exp(dr[index] - variance/2.0); // Exponentiate the field
                    std::poisson_distribution<int> p_dist(density); // Poisson distribution to sample from
                    int num_gals = p_dist(gen); // Poisson sampling
                    totalGals += num_gals;
                    
                    // Output the galaxies in that cell, distributed in a uniform random way within the cell.
                    for (int gal = 0; gal < num_gals; ++gal) {
                        double x = x_min + dist(gen)*delr.x;
                        double y = y_min + dist(gen)*delr.y;
                        double z = z_min + dist(gen)*delr.z;
                        vec3<double> cart = {x + r_min.x, y + r_min.y, z + r_min.z};
                        vec3<double> equa = get_equatorial(cart, cosmo);
                        galaxies[0].push_back(equa.x);
                        galaxies[1].push_back(equa.y);
                        galaxies[2].push_back(equa.z);
                        galaxies[3].push_back(nbar);
                        galaxies[4].push_back(1.0);
                        galaxies[5].push_back(1.0);
                        galaxies[6].push_back(1.0);
                        galaxies[7].push_back(1.0);
                    }
                }
            }
        }
        write_fits_file(out_file, galaxies, col_names, col_forms, col_units);
    }
    
    if (type == density_field) {
        std::vector<double> delta(N.x*N.y*N.z);
        vec3<double> pk_nbw = {0.0, 0.0, 0.0};
        vec3<double> bk_nbw = {0.0, 0.0, 0.0};
        // Loop over the grid and Poisson sample to create the catalog
        //     fout.open(out_file.c_str(), std::ios::app); // Open output file in append mode for multi-tracer case
        //     fout.precision(std::numeric_limits<double>::digits10); // Output full precision
        for (int i = 0; i < N.x; ++i) {
            double x_min = i*delr.x;
            for (int j = 0; j < N.y; ++j) {
                double y_min = j*delr.y;
                for (int k = 0; k < N.z; ++k) {
                    double z_min = k*delr.z;
                    size_t index = k + N.z*(j + N.x*i);
                    
                    double density = n*exp(dr[index] - variance/2.0); // Exponentiate the field
                    std::poisson_distribution<int> p_dist(density); // Poisson distribution to sample from
                    int num_gals = p_dist(gen); // Poisson sampling
                    totalGals += num_gals;
                    // NOTE: Currently this assumes equal weights with those weights set to 1.
                    pk_nbw.x += num_gals;
                    pk_nbw.y += num_gals;
                    pk_nbw.z += nbar*num_gals;
                    
                    bk_nbw.x += num_gals;
                    bk_nbw.y += nbar*num_gals;
                    bk_nbw.z += nbar*nbar*num_gals;
                    delta[index] = num_gals;
                }
            }
        }
        std::cout << pk_nbw.x << " " << pk_nbw.y << " " << pk_nbw.z << "\n";
        std::cout << bk_nbw.x << " " << bk_nbw.y << " " << bk_nbw.z << std::endl;
        write_density_file(out_file, delta, N, L, r_min, pk_nbw, bk_nbw);
    }
    
    if (type == lnknlog) {
        std::ofstream fout(out_file);
        fout.precision(std::numeric_limits<double>::digits10);
        for (int i = 0; i < N.x; ++i) {
            double x_min = i*delr.x;
            for (int j = 0; j < N.y; ++j) {
                double y_min = j*delr.y;
                for (int k = 0; k < N.z; ++k) {
                    double z_min = k*delr.z;
                    size_t index = k + N.z*(j + N.x*i);
                    
                    double density = n*exp(dr[index] - variance/2.0); // Exponentiate the field
                    std::poisson_distribution<int> p_dist(density); // Poisson distribution to sample from
                    int num_gals = p_dist(gen); // Poisson sampling
                    
                    for (int gal = 0; gal < num_gals; ++gal) {
                        double x = x_min + dist(gen)*delr.x;
                        double y = y_min + dist(gen)*delr.y;
                        double z = z_min + dist(gen)*delr.z;
                        fout << x << " " << y << " " << z << " " << nbar << " " << b << " " << 1.0 << "\n";
                        totalGals++;
                    }
                }
            }
        }
        fout.close();
    }
    
    std::cout << "Total number of galaxies: " << totalGals << std::endl;
}
