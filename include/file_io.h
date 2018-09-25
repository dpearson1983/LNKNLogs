#ifndef _FILE_IO_H_
#define _FILE_IO_H_

#include <vector>
#include <string>
#include "tpods.h"

enum FileType{
    dr12,
    patchy,
    density_field,
    lnknlog
};

std::string filename(std::string base, int digits, int num, std::string ext);

FileType setFileType(std::string typeString);

void get_DR12_column_info(std::vector<std::string> &col_names, std::vector<std::string> &col_forms, 
                          std::vector<std::string> &col_units);

void get_DR12_rans_column_info(std::vector<std::string> &col_names, std::vector<std::string> &col_forms, 
                          std::vector<std::string> &col_units);

int write_fits_file(std::string file, std::vector<std::vector<double>> &columns, std::vector<std::string> &col_names,
                     std::vector<std::string> &col_forms, std::vector<std::string> &col_units);

void write_density_file(std::string file, std::vector<double> &delta, vec3<int> N, vec3<double> L, 
                        vec3<double> r_min, vec3<double> pk_nbw, vec3<double> bk_nbw);

#endif
