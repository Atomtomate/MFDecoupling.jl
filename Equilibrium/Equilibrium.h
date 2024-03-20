//
// Created by irakli on 11.01.2023.
//

#include <iostream>
#include <fstream>
#include <complex>
#include <iomanip>
#include </scratch/projects/hhp00048/MFDecoupling/FullGridTE/Equilibrium/eigen-3.4.0/Eigen/Eigenvalues>

using Eigen::MatrixXd;
using Eigen::MatrixXf;
using Eigen::VectorXd;

#ifndef SIAM_TIME_EVOLUTION_EQUILIBRIUM_H
#define SIAM_TIME_EVOLUTION_EQUILIBRIUM_H


void Equilibrium_Selfconsistency(int N_max_eq, double eta, int L, int r, int Sym, double T, double U, double *mu, double ** J, double ** eps, double *eps_d, double *fsigma0, double *gsigma0, double ** rho, MatrixXd& CK_up, MatrixXd& CK_do);

void Impurity(const double U, double const * const mu, double const * const eps, const double T,  double const * const gsigma0, double** rho);
void Chain(const int L, const int r, const int Sym, const double mu, double const * const J, double const * const eps, const double f0, MatrixXd& CK);
void FFF(const int  LL, const double * const EE, const double T, double * const FF);

// For small systems! For Test
void ED_Chain(int L, int r, int Sym,  double mu, const double * J, const double * eps, double f0, MatrixXd& CK);
void Decima_to_binary(double Vect, Eigen::VectorXi& State,int LL);
double Binary_to_decimal(Eigen::VectorXi& State, int LL);

#endif //SIAM_TIME_EVOLUTION_EQUILIBRIUM_H






/*
void Diag(double **matrix, int L, double *eigenvalues,  double *eigenvectors);
extern "C"
{
void dspev_(char* JOBZ, char* UPLO, int *LDZ1, double * AP, double * eigenvalues, double * eigenvectors, int *LDZ2, double *WORK, int *INFO);
}
 */
