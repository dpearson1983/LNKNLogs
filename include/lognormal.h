#ifndef _LOGNORMAL_H_
#define _LOGNORMAL_H_

std::vector<double> fft_freq(int N, double L);

void fill_initial_Pk_grid(std::string in_pk_file, vec3<int> N, vec3<double> L, std::vector<double> &kx, 
                          std::vector<double> &ky, std::vector<double> &kz, fftw_complex *dk, 
                          double b, double f);

std::vector<double> get_dk_initial(vec3<int> N, vec3<double> L, std::string wisdom_file, int nthreads,
                    std::vector<double> &kx, std::vector<double> &ky, std::vector<double> &kz, 
                    double b, double f, std::string in_pk_file);

void get_dk_realization(fftw_complex *dk, fftw_complex *dk_i, 
                        std::vector<double> &kx, std::vector<double> &ky, std::vector<double> &kz,
                        vec3<int> N, vec3<double> L);

void scale_dk_realization(fftw_complex *dk_1, fftw_complex *dk_2, 
                          std::vector<double> &ratio);

void get_galaxies_from_dr(std::vector<double> &dr, vec3<int> N, vec3<double> L, double b, double nbar,
                          std::string out_file);

#endif
