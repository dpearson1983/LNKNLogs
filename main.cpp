/* LNKNLogs v. 1.0
 * David W. Pearson & Lado Samushia
 * Copyright (C) 2016, 2017
 * 
 * This code creates simple lognormal mock catalogs with plane-parallel anisotropies. While these mocks
 * are not ideal for matching surveys, they do produce measurable power spectrum multipoles with 
 * subpercent level accuracy.
 * 
 * GPLv3:
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
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <fftw3.h>
#include <omp.h>
#include "include/tpods.h"
#include "include/harppi.h"
#include "include/lognormal.h"
#include "include/cosmology.h"
#include "include/file_io.h"

int main(int argc, char *argv[]) {
    std::cout << "LNKNLogs v1.0" << std::endl;
    // Check that a parameter file was passed as a command line argument, and then parse it.
    if (argc == 1) {
        std::cout << "Usage: " << std::endl;
        std::cout << "    ./LNKNLogs LNKNLogs.params" << std::endl;
        std::cout << "\nNote that LNKNLogs.params can be any plain text file with any" << std::endl;
        std::cout << "extension that conforms to the HARPPI format. See the included " << std::endl;
        std::cout << "documentation for details." << std::endl;
        return 0;
    }
    parameters p(argv[1]);
    p.print();
    
    std::cout << "Doing some initial setup..." << std::endl;
    cosmology cosmo(p.getd("H_0"), p.getd("Omega_M"), p.getd("Omega_L"));
    vec3<int> N = {p.geti("Nx"), p.geti("Ny"), p.geti("Nz")};
    vec3<double> L = {p.getd("Lx"), p.getd("Ly"), p.getd("Lz")};
    vec3<double> r_min = {p.getd("x_min"), p.getd("y_min"), p.getd("z_min")};
    int N_tot = N.x*N.y*N.z;
    int N_pad = N.x*N.y*2*(N.z/2 + 1);
    int N_rft = N.x*N.y*(N.z/2 + 1);
    
    int num_threads;
    if (p.getb("max_threads")) {
        num_threads = omp_get_max_threads();
    } else if (p.checkParam("num_threads")) {
        num_threads = p.geti("num_threads");
    } else {
        std::cout << "WARNING: Parameter max_threads was set to false and num_threads" << std::endl;
        std::cout << "         was not set. Falling back to single threaded mode." << std::endl;
        num_threads = 1;
    }
    
    std::vector<double> kx = fft_freq(N.x, L.x);
    std::vector<double> ky = fft_freq(N.y, L.y);
    std::vector<double> kz = fft_freq(N.z, L.z);
    
    std::cout << "Getting the field for the random draws..." << std::endl;
    std::vector<fftw_complex> dk_i = get_dk_initial(N, L, p.gets("wisdom_file"), num_threads,
                                              kx, ky, kz, p.getd("b"), p.getd("f"), p.gets("in_pk_file"));
    
    std::cout << "Setting up Fourier transform plan..." << std::endl;
    fftw_init_threads();
    std::vector<fftw_complex> dk(N_rft);
    std::vector<double> dr(N_tot);
    
    FileType type = setFileType(p.gets("fileType"));
    
    fftw_import_wisdom_from_filename(p.gets("wisdom_file").c_str());
    fftw_plan_with_nthreads(num_threads);
    fftw_plan dk2dr = fftw_plan_dft_c2r_3d(N.x, N.y, N.z, dk.data(), dr.data(), FFTW_MEASURE);
    fftw_export_wisdom_to_filename(p.gets("wisdom_file").c_str());
    
    for (int mock = p.geti("start_num"); mock < p.geti("start_num") + p.geti("num_mocks"); ++mock) {
        double start = omp_get_wtime();
        std::string mock_file = filename(p.gets("mock_base"), p.geti("digits"), mock, p.gets("mock_ext"));
        std::cout << "    Creating mock: " << mock_file << std::endl;
        
        get_dk_realization(dk, dk_i, kx, ky, kz, N, L);
        
        fftw_execute(dk2dr);
        
        get_galaxies_from_dr(dr, N, L, r_min, p.getd("b"), p.getd("nbar"), cosmo, mock_file, type);
        std::cout << "    Time: " << omp_get_wtime() - start << " s" << std::endl;
    }
    
    fftw_destroy_plan(dk2dr);
    
    return 0;
}
