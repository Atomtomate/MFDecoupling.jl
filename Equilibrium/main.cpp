#include <iostream>
#include <complex>
#include <cstdlib>


#include "Equilibrium.h"
#include "Read-Write.h"

//#include <stdio.h>
//#include <math.h>
//#include <stdlib.h>
//#include <string.h>
//#include <iomanip>
//#include <stdlib.h>

#define Complex std::complex<double>
using namespace std;


int main(int argc, char **argv) {
    //System Geometry.
    //Number of lattice sites.
    int L=100;
    //Site index to which impurity is coupled 1,...,N.  "0" is reserved for Majorana fermion.
    int r=1;
    //Symetry of the System
    int Sym=1;
    //Temperature.
    double T0=0.;
    double T=0.;
    // Hubbard interaction on the impurity site.
    double U0=4;
    double U=6;
    double V0=0.5;
    double g0_inp, f0_inp;
    int fileID;



    std::cout << "argc == " << argc << std::endl;
    std::cout << "argv == " << argc << std::endl;
    for(int i =0; i<argc; i++)
                printf("%s",argv[i]);
                
    if (argc==7){
        U0=(double) atof(argv[1]);
        V0=(double) atof(argv[2]);
        g0_inp=(double) atof(argv[3]);
        f0_inp=(double) atof(argv[4]);
        fileID=(int) atoi(argv[5]);
        L=(int) atoi(argv[6]);
    }
    else{
        std::cout << "ERROR, argc must be 7, but is " << argc << std::endl;
        exit(1);
    }

    //Equilibrium criteria.
    //Convergence
    double eta=0.0000000001;
    int N_max_eq=100000000;

    //Time evolution.
    //Time step.
    double DeltaT=0.01;
    //Evaluation Time
    double T_eval=10;

    // Information about output, What To Print
    int WTP=2;

    //std::string My_File;			// My_File Name of the my file.
    std::ifstream inFile;			// To open and read from the file.




    cout << "Program started! \n" ;
    cout << "Reading from Input_File" << endl;

    //  Here we open Input_File to read input parameters
    inFile.open("Input_File");
    // Check for error
    if (inFile.fail()){
        cerr << "Error opening Input_File" << endl;
        exit(1);
    }

    //Chemical potential
    double *mu;
    mu = new double [2];


    // Here we start to read input parameters:  Model parameters
    Input_File_1(&r, &Sym, mu, inFile);

    MatrixXd CK_up = MatrixXd::Zero(2*L+2,2*L+2);
    MatrixXd CK_do = MatrixXd::Zero(2*L+2,2*L+2);

    //Declaring matrices
    //Hopping J
    // J[L] is equivalent to Vr.
    double **J0, **J;	//[L+1]
    J0 = new double *[2];
    J = new double * [2];
    for (int s=0; s<2; ++s){
        J0[s] = new double [L+1];
        J[s] = new double [L+1];
    }

    //Onsite energies chain
    //Onsite energies for the chain.
    double **eps0, **eps; //[L]
    eps0 = new double * [2];
    eps = new double * [2];
    for (int s=0; s<2; ++s){
        eps0[s] = new double [L+1];
        eps[s] = new double [L+1];
    }


    // onsite energies impurity
    //Onsite energies for the impurity.
    double *eps_d0, *eps_d; //[L]
    eps_d0 = new double [2];
    eps_d = new double [2];

    //Mean-field parameters g and f
    //$g_\sigma=V_r C_{0r,\sigma}$.
    double *g0; // [2]
    double  *f0;	//[2]
    g0 = new double [2];
    f0 = new double [2];

    //rho-matrix for the impurity
    double **rho0;	//[4][4]
    //Complex **rho;	//[4][4]
    rho0 = new double * [4];
    for(int i=0; i<4; ++i )
        rho0[i]= new double[4];


    // Here we continue to read input parameters
    Input_File_2(Sym, L, mu, &T0, eps_d0, eps0, J0, &eta, &N_max_eq, &WTP,  inFile);

    J0[0][L]=V0;
    J0[1][L]=V0;

    for(int sigma=0; sigma <2; ++sigma){
        eps_d0[sigma]= eps_d0[sigma] + mu[sigma] +U0/2.;
    }

    //  Here we close Input_File
    inFile.close();
    for(int sigma=0; sigma<2; ++sigma){
        g0[sigma]=g0_inp;
        f0[sigma]=f0_inp;
        cout << sigma <<"\t"<<g0[sigma]<<"\t"<<f0[sigma]<<endl;
    }

    //Equilibrium Calculations
    //  In input file r runs from 1 to L, while in Hamiltonian chain coordinates run from 0 to L-1
    //  Therefore in input I wrote r-1 to take into account this shift
    Equilibrium_Selfconsistency(N_max_eq, eta, L, r-1, Sym, T0, U0, mu, J0, eps0, eps_d0, f0, g0, rho0, CK_up, CK_do);
    cout << "Equilibrium calculations done, writing to ID=" << fileID << endl;
    Output_File_equilibrium(L, CK_up, CK_do,  rho0, f0, g0, U0, J0[0][L], WTP, fileID);


    /* Eigen::MatrixXd rho_Xd = Eigen::MatrixXd::Zero(4, 4); */
    /* for (int i = 0; i < 4; ++i) { */
    /*     for (int j = i; j < 4; ++j) { */
    /*         rho_Xd(i, j) = rho0[i][j]; */
    /*         rho_Xd(j, i) = rho0[j][i]; */
    /*     } */
    /* } */
    /* res.rho0 = rho_Xd; */
    /* res.CK_up = CK_up; */
}
