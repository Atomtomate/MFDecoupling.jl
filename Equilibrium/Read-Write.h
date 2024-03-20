//
// Created by Irakli on 30.01.23.
//

#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>
#include </scratch/projects/hhp00048/MFDecoupling/FullGridTE/Equilibrium/eigen-3.4.0/Eigen/Eigenvalues>

using namespace std;

#ifndef SIAM_TIME_EVOLUTION_READ_WRITE_H
#define SIAM_TIME_EVOLUTION_READ_WRITE_H



void Input_File_1(int *r, int *Sym,double *mu, std::ifstream &inFile);
void Input_File_2(int Sym, int L, const double *mu, double *T0, double *eps_d0, double **eps0, double **J0, double *eta, int *N_max_eq, int *WTP, std::ifstream &inFile);

template <typename T>
std::string to_string_with_precision(const T& a_value, const int n = 6)
{
    std::ostringstream out;
    out << std::setprecision(n) << a_value;
    return out.str();
}

void Output_File_equilibrium(int L, Eigen::MatrixXd CK_up, Eigen::MatrixXd CK_do,  double **rho0, double *f0, double *g0, double U0, double V0, int WTP, int);

#endif //SIAM_TIME_EVOLUTION_READ_WRITE_H
