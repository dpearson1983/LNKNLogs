#include <fstream>
#include <CCfits/CCfits>
#include <vector>
#include <string>
#include "../include/tpods.h"
#include "../include/file_io.h"

std::string filename(std::string base, int digits, int num, std::string ext) {
    std::stringstream file;
    file << base << std::setw(digits) << std::setfill('0') << num << "." << ext;
    return file.str();
}

void get_DR12_column_info(std::vector<std::string> &col_names, std::vector<std::string> &col_forms, 
                          std::vector<std::string> &col_units) {
    col_names.push_back("RA");
    col_names.push_back("DEC");
    col_names.push_back("Z");
    col_names.push_back("NZ");
    col_names.push_back("WEIGHT_FKP");
    col_names.push_back("WEIGHT_SYSTOT");
    col_names.push_back("WEIGHT_NOZ");
    col_names.push_back("WEIGHT_CP");
    
    col_forms.push_back("D");
    col_forms.push_back("D");
    col_forms.push_back("E");
    col_forms.push_back("E");
    col_forms.push_back("E");
    col_forms.push_back("E");
    col_forms.push_back("E");
    col_forms.push_back("E");
    
    col_units.push_back("DEG");
    col_units.push_back("DEG");
    col_units.push_back("");
    col_units.push_back("Mpc^-1");
    col_units.push_back("");
    col_units.push_back("");
    col_units.push_back("");
    col_units.push_back("");
}

void get_DR12_rans_column_info(std::vector<std::string> &col_names, std::vector<std::string> &col_forms, 
                          std::vector<std::string> &col_units) {
    col_names.push_back("RA");
    col_names.push_back("DEC");
    col_names.push_back("Z");
    col_names.push_back("NZ");
    col_names.push_back("WEIGHT_FKP");
    
    col_forms.push_back("D");
    col_forms.push_back("D");
    col_forms.push_back("E");
    col_forms.push_back("E");
    col_forms.push_back("E");
    
    col_units.push_back("DEG");
    col_units.push_back("DEG");
    col_units.push_back("");
    col_units.push_back("Mpc^-1");
    col_units.push_back("");
}

FileType setFileType(std::string typeString) {
    FileType type;
    if (typeString == "DR12") {
        type = dr12;
    } else if (typeString == "Patchy") {
        type = patchy;
    } else if (typeString == "density_field") {
        type = density_field;
    } else if (typeString == "lnknlog") {
        type = lnknlog;
    }
    return type;
}

int write_fits_file(std::string file, std::vector<std::vector<double>> &columns, std::vector<std::string> &col_names,
                     std::vector<std::string> &col_forms, std::vector<std::string> &col_units) {
    std::unique_ptr<CCfits::FITS> pFits(new CCfits::FITS(file, CCfits::Write));
    
    std::string hduName("LNKNLOG");
    
    unsigned long rows = columns[0].size();
    
    CCfits::Table *newTable = pFits->addTable(hduName, rows, col_names, col_forms, col_units);
    for (int i = 0; i < col_names.size(); ++i) {
        newTable->column(col_names[i]).write(columns[i],1);
    }
    
    return 1;
}

void write_density_file(std::string file, std::vector<double> &delta, vec3<int> N, vec3<double> L, 
                        vec3<double> r_min, vec3<double> pk_nbw, vec3<double> bk_nbw) {
    std::ofstream fout(file, std::ios::out|std::ios::binary);
    fout.write((char *) &L, 3*sizeof(double));
    fout.write((char *) &r_min, 3*sizeof(double));
    fout.write((char *) &pk_nbw, 3*sizeof(double));
    fout.write((char *) &bk_nbw, 3*sizeof(double));
    fout.write((char *) delta.data(), N.x*N.y*N.z*sizeof(double));
    fout.close();
}
