//
// Created by Irakli on 30.01.23.


#include "Read-Write.h"
using namespace std;

//Model parameters
//=====================================================================================================================================================================================
void Input_File_1(int *r, int *Sym, double *mu, std::ifstream &inFile){
    std::string dummy_text;			// dummy_text - To save dummy information.
    double mu_inp, h_inp;
    int L_inp, r_inp, Sym_inp;

    // Skipping line: // Model parameters
    getline(inFile,dummy_text);

    // Skipping line: //  L		*integer*	Number of sites in chain.
    getline(inFile,dummy_text);
    // Reading L
    inFile >> L_inp;
    //*L=L_inp;
    if( inFile.fail() && !inFile.eof() ){
        //clear error state
        inFile.clear();
        //set x to NaN
        cout << "Set up for parameters is incorrect for L" << endl;
        exit(1);
    }
    getline(inFile,dummy_text); // moving to the next line

    // Skipping line: //  r		*integer*	Site coordinate to which impurity is coupled. 0, ... , L-1.
    getline(inFile,dummy_text);
    // Reading r
    inFile >> r_inp;
    *r=r_inp;
    if( inFile.fail() && !inFile.eof() ){
        //clear error state
        inFile.clear();
        //set x to NaN
        cout << "Set up for parameters is incorrect for r" << endl;
        exit(1);
    }
    getline(inFile,dummy_text); // moving to the next line

    // Skipping line: //  Sym		*integer*	Symmetry of the system:
    getline(inFile,dummy_text);
    // Skipping line: //				(0) No symmetry; (1) All $J$ are equal OBC; (2) All $J$ are equal PBC; (3) Ionic chain OBC; (4) Ionic chain PBC.
    getline(inFile,dummy_text);
    // Reading Sym
    inFile >> Sym_inp;
    *Sym=Sym_inp;
    if( inFile.fail() && !inFile.eof() ){
        //clear error state
        inFile.clear();
        //set x to NaN
        cout << "Set up for parameters is incorrect for Sym" << endl;
        exit(1);
    }
    if ((*Sym ==3 || *Sym ==4) ){ //&& *L % 2 !=0 
        cout << "Wrong output! For ionic chains the number of sites should be even. Calculations are aborted!"<<endl;
       exit(1);
    }
    getline(inFile,dummy_text); // moving to the next line

    // Skipping line: // mu		*double*	Chemical potential.  $mu_up=mu-h/2$ and $mu_down=mu+h/2$ .
    getline(inFile,dummy_text);
    // Reading mu
    inFile >> mu_inp;
    if( inFile.fail() && !inFile.eof() ){
        //clear error state
        inFile.clear();
        //set x to NaN
        cout << "Set up for parameters is incorrect for mu" << endl;
        exit(1);
    }
    getline(inFile,dummy_text); // moving to the next line


    // Skipping line: // h		*double*	Magnetic field.
    getline(inFile,dummy_text);
    // Reading h
    inFile >> h_inp;
    if( inFile.fail() && !inFile.eof() ){
        //clear error state
        inFile.clear();
        //set x to NaN
        cout << "Set up for parameters is incorrect for h" << endl;
        exit(1);
    }
    getline(inFile,dummy_text); // moving to the next line

    mu[0]=mu_inp-0.5*h_inp;
    mu[1]=mu_inp+0.5*h_inp;

}


//=====================================================================================================================================================================================
void Input_File_2(const int Sym, const int L, const double *mu, double *T0, double *eps_d0, double **eps0, double **J0, double *eta, int *N_max_eq, int *WTP, std::ifstream &inFile){
    std::string dummy_text;			// dummy_text - To save dummy information.
    double epsd0_up, epsd0_down, epsd_up, epsd_down;
    double eps0_1, eps0_2, eps_1, eps_2;
    double J0_1, J0_2, V0_r, J_1, J_2, V_r;
    double g0_inp, f0_inp;
    double eta_inp;
    int N_max_eq_inp;
    double T0_inp, U0_inp, T_inp, U_inp;
    double DeltaT_inp, T_eval_inp;
    int WTP_inp;

    std::ifstream inFile_0;			// To open and read from the file.



    // Skipping line: // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    getline(inFile,dummy_text);
    // Skipping line: // Before the quench !!!!!
    getline(inFile,dummy_text);
    // Skipping line: // Physical parameters:
    getline(inFile,dummy_text);

    // Skipping line: // T0		*double*	Temperature of the system.
    getline(inFile,dummy_text);
    // Reading T0
    inFile >> T0_inp;
    *T0=T0_inp;
    if( inFile.fail() && !inFile.eof() ){
        //clear error state
        inFile.clear();
        //set x to NaN
        cout << "Set up for parameters is incorrect for T0" << endl;
        exit(1);
    }
    getline(inFile,dummy_text); // moving to the next line

    /*
    // Skipping line: // U0		*double*	Hubbard interaction on the impurity site.
    getline(inFile,dummy_text);
    // Reading U0
    inFile >> U0_inp;
    *U0=U0_inp;
    if( inFile.fail() && !inFile.eof() ){
        //clear error state
        inFile.clear();
        //set x to NaN
        cout << "Set up for parameters is incorrect for U0" << endl;
        exit(1);
    }
    getline(inFile,dummy_text); // moving to the next line
    */

    // Skipping line: // eps_d0		*double*    Onsite energy on the impurity site. !!! The first number is the correction for spin up and the second for spin down !!!
    getline(inFile,dummy_text);
    // Skipping line: //			$\varepsilon_{d,\sigma}=U/2+mu_\sigma - \varepsilon_{input,\sigma}
    getline(inFile,dummy_text);
    // Skipping line: //				!!! For half-filled impurity we should have 0 0 !!!
    getline(inFile,dummy_text);
    // Reading eps_d0
    inFile >> epsd0_up;
    inFile >> epsd0_down;
    if( inFile.fail() && !inFile.eof() ){
        //clear error state
        inFile.clear();
        //set x to NaN
        cout << "Set up for parameters is incorrect for eps_d0" << endl;
        exit(1);
    }
    getline(inFile,dummy_text); // moving to the next line

    eps_d0[0]=epsd0_up;   //It will be corrected later
    eps_d0[1]=epsd0_down;   //It will be corrected later


    // Skipping line: // eps0		*double*	Onsite energy for the chain.
    getline(inFile,dummy_text);
    // Skipping line: //				If Sym (0) program will read from separate file; For (1) and (2) it will be automatically zeo! For (3) and (4) will be two numbers (for the first and the second site of the unit cell).
    getline(inFile,dummy_text);
    // Reading eps0
    inFile >> eps0_1;
    inFile >> eps0_2;
    if( inFile.fail() && !inFile.eof() ){
        //clear error state
        inFile.clear();
        //set x to NaN
        cout << "Set up for parameters is incorrect for eps0" << endl;
        exit(1);
    }
    getline(inFile,dummy_text); // moving to the next line


    if (Sym==0){
        inFile_0.open("Input_eps0");
         /*
        // Check for error
        if (inFile_0.fail()){
            cerr << "Error opening Input_eps0" << endl;
            exit(1);
        }
        */
         int read_int;
         for (int i=0; i<L; ++i){
             inFile_0 >> read_int;
             inFile_0 >> eps0[0][i];
             eps0[1][i] = eps0[0][i];
             /*
             if( inFile_0.fail() && !inFile_0.eof() ){
                 //clear error state
                 inFile_0.clear();
                 //set x to NaN
                 cout << "Set up for parameters is incorrect for eps0 in Input_eps0!" << endl;
                 exit(1);
             }
             */
             getline(inFile_0,dummy_text); // moving to the next line
         }
        inFile_0.close();
    }
    else if (Sym ==1 || Sym ==2){
        for (int i=0; i<L; ++i)
            for(int sigma=0; sigma <2; ++sigma)
                eps0[sigma][i]=0.0;
    }
    else if (Sym ==3 || Sym ==4){
        for (int i=0; i<L/2; ++i)
            for(int sigma=0; sigma <2; ++sigma){
                eps0[sigma][2*i]=eps0_1;
                eps0[sigma][2*i+1]=eps0_2;
            }
    }

    // Skipping line: // J0		*double*	Hopping amplitude between neighboring sites.
    getline(inFile,dummy_text);
    // Skipping line: //				If Sym (0) program will read from separate file; For (1) and (2) there will be one number; For (3) and (4) will be two numbers (inside the unit cell and between the unit cells).
    getline(inFile,dummy_text);
    // Reading J0
    inFile >> J0_1;
    if (Sym == 3 || Sym == 4 ){
        inFile >> J0_2;
    }
    if( inFile.fail() && !inFile.eof() ){
        //clear error state
        inFile.clear();
        //set x to NaN
        cout << "Set up for parameters is incorrect for J0" << endl;
        exit(1);
    }
    getline(inFile,dummy_text); // moving to the next line

    if (Sym==0){
        inFile_0.open("Input_J0");
         /*
        // Check for error
        if (inFile_0.fail()){
            cerr << "Error opening Input_J0" << endl;
            exit(1);
        }
        */
         int read_int;
         for (int i=0; i<L; ++i){
             inFile_0 >> read_int;
             inFile_0 >> J0[0][i];
             J0[1][i] = J0[0][i];
             /*
             if( inFile_0.fail() && !inFile_0.eof() ){
                 //clear error state
                 inFile_0.clear();
                 //set x to NaN
                 cout << "Set up for parameters is incorrect for J0 in Input_J0!" << endl;
                 exit(1);
             }
             */
             getline(inFile_0,dummy_text); // moving to the next line
         }
        inFile_0.close();

    }
    else if (Sym ==1){
        for (int sigma = 0; sigma < 2; ++sigma){
            for (int i=0; i<L-1; ++i)
                J0[sigma][i] = J0_1;
            J0[sigma][L-1] = 0.;
        }
    }
    else if (Sym ==2){
        for (int i=0; i<L; ++i)
            for(int sigma=0; sigma <2; ++sigma)
                J0[sigma][i]=J0_1;
    }
    else if (Sym ==3){
        for (int i=0; i<L/2; ++i)
            for(int sigma=0; sigma <2; ++sigma){
                J0[sigma][2*i]=J0_1;
                J0[sigma][2*i+1]=J0_2;
            }
        for(int sigma=0; sigma <2; ++sigma)
            J0[sigma][L-1]=0.;
    }
    else if (Sym ==4){
        for (int i=0; i<L/2; ++i)
            for(int sigma=0; sigma <2; ++sigma){
                J0[sigma][2*i]=J0_1;
                J0[sigma][2*i+1]=J0_2;
            }
    }

    /*
    // Skipping line: // Vr0		*double*	Coupling between impurity and the chain.
    getline(inFile,dummy_text);
    // Skipping line: //				If Sym (0)-(4) there will be one number.
    getline(inFile,dummy_text);
    // Reading Vr0
    inFile >> V0_r;
    if( inFile.fail() && !inFile.eof() ){
        //clear error state
        inFile.clear();
        //set x to NaN
        cout << "Set up for parameters is incorrect for Vr0" << endl;
        exit(1);
    }
    getline(inFile,dummy_text); // moving to the next line
    */

    for(int sigma=0; sigma <2; ++sigma)
        J0[sigma][L]=0.5; //V0_r;    //It will be corrected later


    // Skipping line: // Convergence parameters
    getline(inFile,dummy_text);


    // Skipping line: // eta		*double*	Convergence accuracy in the equilibrium state
    getline(inFile,dummy_text);
    // Reading eta
    inFile >> eta_inp;
    *eta=eta_inp;
    if( inFile.fail() && !inFile.eof() ){
        //clear error state
        inFile.clear();
        //set x to NaN
        cout << "Set up for parameters is incorrect for eta" << endl;
        exit(1);
    }
    getline(inFile,dummy_text); // moving to the next line

    // Skipping line: // N_max_eq	*double*	Maximum iterations in the convergence loop at equilibrium
    getline(inFile,dummy_text);
    // Reading N_max_eq
    inFile >> N_max_eq_inp;
    *N_max_eq=N_max_eq_inp;
    if( inFile.fail() && !inFile.eof() ){
        //clear error state
        inFile.clear();
        //set x to NaN
        cout << "Set up for parameters is incorrect for eta" << endl;
        exit(1);
    }
    getline(inFile,dummy_text); // moving to the next line



    // Skipping line: // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    getline(inFile,dummy_text);

    // Skipping line: // Information about output
    // Skipping line: //WTP		*integer*	 (0) Print everything.  I might add other options later!!!
    getline(inFile,dummy_text);
    getline(inFile,dummy_text);
    // Reading WTP
    inFile >> WTP_inp;
    *WTP = WTP_inp;
    if( inFile.fail() && !inFile.eof() ){
        //clear error state
        inFile.clear();
        //set x to NaN
        cout << "Set up for parameters is incorrect for WTP" << endl;
        exit(1);
    }
    getline(inFile,dummy_text); // moving to the next line



}


//1234567890
//===============================================================================================
void Output_File_equilibrium(int L, Eigen::MatrixXd CK_up, Eigen::MatrixXd CK_do,  double **rho0, double *f0, double *g0, double U0, double V0,int WTP, int fileID) {

    std::ofstream outFile;            // To open and write in the file.

    CK_do = CK_up;
    double Nup_c = 0;
    double Ndown_c = 0;
    for (int i = 0; i < L; ++i) {
        Nup_c += CK_up(i, i);
        Ndown_c += CK_do(i, i);
    }
    Nup_c = Nup_c / L;
    Ndown_c = Ndown_c / L;


    cout << "Filling of the chain is Nup_c=" << Nup_c << "\t Ndown_c=" << Ndown_c << endl;

    double Nup_imp = rho0[0][0] + rho0[1][1];
    double Ndown_imp = rho0[0][0] + rho0[2][2];


    cout << "Filling of the impurity is Nup_imp=" << Nup_imp << "\t Ndown_imp=" << Ndown_imp << endl;


    /* outFile.open("Main_Results_U" + to_string_with_precision(U0, 5) + "V" + to_string_with_precision(V0, 5) + ".dat", ios::out); */
    /* outFile << g0[0] << "\t" << f0[0] << "\t"; */
    /* for (int i = 0; i < 4; ++i) { */
    /*     for (int j = i; j < 4; ++j) { */
    /*         outFile << std::setprecision (17) << rho0[i][j] << "\t"; */
    /*     } */
    /* } */
    /* outFile << endl; */
    /* outFile.close(); */


    if (WTP > 0) {
        outFile.open("CK_" + std::to_string(fileID) + ".dat", ios::out);
        outFile << CK_up << endl;
        outFile.close();

        /*
         outFile.open("CK_up.dat");
         outFile << CK_up<<endl;
         outFile.close();

         outFile.open("CK_do.dat");
         outFile << std::setprecision (17) << CK_up<<endl;
         outFile.close();
         */


        Eigen::MatrixXd rho_Xd = Eigen::MatrixXd::Zero(4, 4);

        for (int i = 0; i < 4; ++i) {
            for (int j = i; j < 4; ++j) {
                rho_Xd(i, j) = rho0[i][j];
                rho_Xd(j, i) = rho0[j][i];
            }
        }


        outFile.open("rho0_" + std::to_string(fileID) + ".dat", ios::out);
        outFile << std::setprecision (17) << rho_Xd << endl;
        outFile.close();
    }
}
