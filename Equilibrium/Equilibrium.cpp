//
// Created by Irakli on 11.01.2023.
//

#include "Equilibrium.h"

#define Complex std::complex<double>

using namespace std;

void Equilibrium_Selfconsistency(int N_max_eq, double eta, int L, int r, int Sym, double T, double U, double *mu, double ** J, double ** eps, double *eps_d, double *fsigma0, double *gsigma0, double ** rho, MatrixXd& CK_up, MatrixXd& CK_do){

    double diff=10.;
    int count=0;
    double *ftemp, *gtemp;

    ftemp = new double [2];
    gtemp = new double [2];
    cout << "Equilibrium calculations are started"<<endl;

    while (diff > eta) {
        count += 1;

        //Impurity
        Impurity(U, mu, eps_d, T, gsigma0, rho);

        // Chain: spin up fermions
        Chain(L, r, Sym, mu[0], J[0], eps[0],fsigma0[0],CK_up);
        //cout << CK_up <<endl;
        // Chain: spin down fermions
        //Chain(L, r, Sym, mu[1], J[1], eps[1],fsigma0[1],CK_do);   !!!!!
        //cout << CK_do <<endl;

        for (int sigma = 0; sigma < 2; ++sigma) {
            ftemp[sigma] = fsigma0[sigma];
            gtemp[sigma] = gsigma0[sigma];
        }


        // Matrice indexes are shifted by one in contrast to PDF !!!
        fsigma0[0] = J[0][L] * (rho[0][2] + rho[1][3]);
        fsigma0[1] = J[1][L] * (rho[0][1] + rho[2][3]);

        gsigma0[0] = J[0][L] * ( CK_up(0,2*L) + CK_up(0,2*L+1));
        gsigma0[1] = gsigma0[0]; // J[1][L] * ( CK_do(0,2*L) + CK_do(0,2*L+1)); !!!!!

        diff = 0.0;
        for (int sigma = 0; sigma < 2; ++sigma) {
            diff += (fabs(ftemp[sigma]) - fabs(fsigma0[sigma])) * (fabs(ftemp[sigma]) - fabs(fsigma0[sigma]));
            diff += (fabs(gtemp[sigma]) - fabs(gsigma0[sigma])) * (fabs(gtemp[sigma]) - fabs(gsigma0[sigma]));
        }
        diff = sqrt(diff / 4.);
        cout << "After " << count << " iteration diff=" << diff << "\t" << "fsigma0[0]=" << fsigma0[0] << " fsigma0[1]=" << fsigma0[1] << "\t" << "gsigma0[0]=" << gsigma0[0] << " gsigma0[1]=" << gsigma0[1] << endl;

        if (count > N_max_eq) {
            cout << "Convergence after " << N_max_eq << "iterations do not reach desired accuracy eta=" << eta << ".\t  We obtain diff=" << diff << endl;
            return;
        }

    }

    CK_do=CK_do;

    //    if (count<N_max_eq+1) {
    cout << "Convergence is reached after " << count << " iterations. We obtain diff=" << diff << endl;
    //  }
    cout << "RESULT: f = " << fsigma0[1] << " // g = " << gsigma0[1] << endl;

    cout << "Equilibrium calculations are finished"<<endl;

    free(ftemp);
    free(gtemp);
}


void Chain(const int L, const int r, const int Sym, const double mu, double const * const J, double const * const eps, const double f0, MatrixXd& CK){
    //Hamiltonian
    MatrixXd H = MatrixXd::Zero(2*L+2,2*L+2);
    //Diagonal Matrix
    MatrixXd D = MatrixXd::Zero(2*L+2,2*L+2);
    //Unitary transformation
    MatrixXd U = MatrixXd::Zero(2*L+2,2*L+2);
    // Eigenvalues
    VectorXd El= VectorXd::Zero(2*L+2);

    //MatrixXd CK = MatrixXd::Zero(2*L+2,2*L+2);

    //Diagonal Terms
    for(int i=0; i<L; ++i){
        H(i,i)= eps[i] - mu;
        H(i+L,i+L)=-H(i,i);
    }
    //Off-diagonal Terms
    H(0,L-1)=J[L-1];
    H(L-1,0)=J[L-1];
    H(L,L-1+L)=-J[L-1];
    H(L-1+L,L)=-J[L-1];
    for(int i=0; i<L-1; ++i){
        H(i,i+1)= J[i];
        H(i+1,i)= J[i];
        H(i+L,i+1+L)=-J[i];
        H(i+1+L,i+L)=-J[i];
    }


    //coupling to Majorana
    H(r,2*L)=f0;
    H(r,2*L+1)=f0;
    H(2*L,r)=f0;
    H(2*L+1,r)=f0;
    H(r+L,2*L)=-f0;
    H(r+L,2*L+1)=-f0;
    H(2*L,r+L)=-f0;
    H(2*L+1,r+L)=-f0;

    //cout << H << endl;



    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(H);

    El= es.eigenvalues();
    U = es.eigenvectors();

    for (int i=0; i<2*L+2; ++i){
        if (El(i)<-1e-10)
            D(i,i)=1.;
        else if (El(i)>1e-10)
            D(i,i)=0.;
        else
            D(i,i)=0.5;
    }


/*
    MatrixXd M1 = MatrixXf::Random(2*L+2,2*L+2);
    cout << M1 <<endl;
    Eigen::Map<MatrixXf,0,Eigen::OuterStride<> > M2(M1.data(), M1.rows()-L, M1.cols(), Eigen::OuterStride<>(M1.outerStride()*1));
    cout << "1 column over 3:" << endl << M2 << "\n";
    cout << endl;

    cout << M2*M1;
*/

    CK =U*D* U.transpose();

}



void Impurity(const double U, double const * const mu, double const * const eps, const double T,  double const * const gsigma0, double** rho){
    double **H;
    // Two-dimensional matrix is written as one dimensional
    double **UU;
    double *El;
    double *FF;

    //Declaring matrices
    //Hamiltonian
    H = new double *[4];
    for (int i = 0; i < 4; ++i)
        H[i] = new double[4];

    El  = new double[4];
    UU= new double * [4];
    for (int i=0; i<4; ++i)
        UU[i] = new double [4];

    FF = new double [4];

    // Filling Hamiltonian
    H[0][0]=U+mu[0]+mu[1]-eps[0]-eps[1];
    H[1][1]=mu[0]-eps[0];
    H[2][2]=mu[1]-eps[1];
    H[3][3]=0;

    H[0][1]=gsigma0[1];
    H[1][0]=gsigma0[1];
    H[0][2]=gsigma0[0];
    H[2][0]=gsigma0[0];
    H[0][3]=0.;
    H[3][0]=0.;
    H[1][2]=0.;
    H[2][1]=0.;
    H[1][3]=gsigma0[0];
    H[3][1]=gsigma0[0];
    H[2][3]=gsigma0[1];
    H[3][2]=gsigma0[1];
    
    //Finding Eigenvalues and eigenvectors
    Eigen::MatrixXd A (4, 4);
    for (int i=0; i<4; ++i){
        for (int j=0; j<4; ++j){
            A(i,j)=H[i][j];
        }
    }
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(A);
    const auto es_tmp = es.eigenvalues();
    for (int i=0; i<4; ++i){
        El[i]=es_tmp[i];
    }

    const auto ev_tmp = es.eigenvectors();
    for (int i=0; i<4; ++i)
        for (int j=0; j<4; ++j)
            UU[i][j]=ev_tmp(i,j);

    FFF(4, El, T, FF);

    //In TeX file we have expression $\rho_{ij}=\frac{1}{Z}\sum_{l=1}^{4} \exp(-\frac{E_l}{k_B T}) \alpha_{lj}^*\alpha_{li}
    //Note that \hat \alpha= \hat U^{-1}=U^\dagger$ Here $U$ is transformation matrix diagonalizing the Hamiltonian. In our case $UU$
    //Note by construction (all entries are real numbers) eigenvectors are real numbers!!! So \hat\alpha = U^T
    // \alpha_{li} = U_{il} = eigenvectors[l*4+i]   and   \alpha_{lj}^*=U_{jl}
    for(int i=0; i<4; ++i)
        for(int j=i; j<4; ++j){
            double temp=0.;
            for (int l=0;l <4; ++l) {
                temp += FF[l] * UU[i][l] * UU[j][l];
            }

            rho[i][j] =temp;
            rho[j][i] =temp;
        }
    for(int i=0; i<4; ++i){
        free(UU[i]);
        free(H[i]);
    }
    free(UU);
    free(H);

    free(El);
    free(FF);
}



void FFF(const int  LL, const double * const EE, const double T, double * const FF){

    if (T<0.00000000001) {
        int count = 1;
        for (int l = 1; l < LL; ++l) {
            if (EE[l] - EE[0] < 1e-12) {
                //cout << l << "\t" <<count<<"\t"<< EE[l] - EE[0] << endl;
                count += 1;
            }
            else
                break;
        }

        for (int l = 0; l < count; ++l)
            FF[l] = 1. / count;
        for (int l = count ; l < LL; ++l)
            FF[l] = 0.0;
    }
    else{
        double Z=0.0;
        for (int l=0; l<LL;++l)
            Z+=exp(-1.0*EE[l]/T);

        for (int l=0; l<LL;++l)
            FF[l]=exp(-1.0*EE[l]/T)/Z;
    }

}


// For small systems! For Test
void ED_Chain(int L, int r, int Sym,  double mu, const double * J, const double * eps, double f0, MatrixXd& CK){

    int Lstate=  pow(2,L);
    //Hamiltonian
    MatrixXd H = MatrixXd::Zero(Lstate,Lstate);

    VectorXd El= VectorXd::Zero(Lstate);
    //Diagonal Matrix
    MatrixXd D = MatrixXd::Zero(Lstate,Lstate);
    //Unitary transformation
    MatrixXd U = MatrixXd::Zero(Lstate,Lstate);
    // Eigenvalues
    Eigen::VectorXi Vect = Eigen::VectorXi::Zero(Lstate);
    Eigen::VectorXi State = Eigen::VectorXi::Zero(L+1);



    int istate=1;
    for(int ii=1; ii<2*Lstate-1; ++ii){
        int n=ii;
        int count=0;
        for(int i=0; n>0; i++){
            count += n%2 ;
            n= n/2;
        }
        if (count % 2 == 0){
            Vect(istate)=ii;
            istate +=1;
        }
    }



    // Diagonal elements
    for (int i=0; i<L; ++i){
        for(int ir=0; ir<Lstate; ++ir){
            Decima_to_binary(Vect(ir), State,L+1);
            if (State(i)==1){
                H(ir,ir)+=eps[i] - mu;;
            }
        }
    }


    // Hopping between chain sites
    for (int i=0; i<L-1; ++i){
        for(int ir=0; ir<Lstate; ++ir){
            Decima_to_binary(Vect(ir), State,L+1);
            if (State(i)==0 &&  State(i+1)==1){
                State(i)=1;
                State(i+1)=0;
                double Vector=Binary_to_decimal(State, L+1);

                for(int il=0; il<Lstate;++il){
                    if (Vector==Vect(il)){
                        H(il,ir)+=J[i];
                        H(ir,il)+=J[i];
                        break;
                    }
                }
            }
        }
    }

    // Hopping between chain and Majorana
    for(int ir=0; ir<Lstate; ++ir){
        Decima_to_binary(Vect(ir), State,L+1);
        if (State(0)==1 &&  State(L)==0){

            State(0)=0;
            //Checking sign
            int count=1;
            for(int ii=0; ii<L; ++ii){
                if (State(ii)==1) {
                    count *= -1;
                }
            }
            State(L)=1;

            double Vector=Binary_to_decimal(State, L+1);

            for(int il=0; il<Lstate;++il){
                if (Vector==Vect(il)){
                    H(il,ir)+=count*f0;
                    H(ir,il)+=count*f0;
                    break;
                }
            }
        }
    }

    for(int ir=0; ir<Lstate; ++ir){
        Decima_to_binary(Vect(ir), State,L+1);
        if (State(0)==1 &&  State(L)==1){

            State(0)=0;
            //Checking sign
            int count=1;
            for(int ii=1; ii<L; ++ii){
                if (State(ii)==1) {
                    count *= -1;
                }
            }
            State(L)=0;

            double Vector=Binary_to_decimal(State, L+1);

            for(int il=0; il<Lstate;++il){
                if (Vector==Vect(il)){
                    H(il,ir)+=-count*f0;
                    H(ir,il)+=-count*f0;
                    break;
                }
            }
        }
    }


    //cout << H << endl;


    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(H);

    El= es.eigenvalues();
    U = es.eigenvectors();

    //cout << El << endl;

    //cout <<U<<endl;

    //Calculating correlators C(i,j) i \neq j not Majorana
    for (int i=0; i<L;++i){
        for(int j=i+1; j<L;++j){

            CK(i,j)=0;

            for(int ir=0; ir<Lstate; ++ir){
                Decima_to_binary(Vect(ir), State,L+1);
                int count = 1;
                if (State(j)==1 &&  State(i)==0){


                    //Checking sign
                    for(int ii=i; ii<j; ++ii){
                        if (State(ii)==1) {
                            count *= -1;
                        }
                    }
                    State(j)=0;
                    State(i)=1;

                    double Vector=Binary_to_decimal(State, L+1);

                    // Calculating CK
                    for(int il=0; il<Lstate;++il){
                        if (Vector==Vect(il)){

                            CK(i,j)+=count*U(il,0)*U(ir,0);
                            break;
                        }
                    }
                }
            }

            CK(j,i)=CK(i,j);

            CK(i+L,j+L)=-CK(i,j);
            CK(j+L,i+L)=CK(i+L,j+L);
        }
    }

     //Calculating correlators C(i,i) not Majorana
    for (int i=0; i<L;++i){

        CK(i,i)=0;
        for(int ir=0; ir<Lstate; ++ir) {
            Decima_to_binary(Vect(ir), State, L + 1);
            if (State(i) == 1) {
                // Calculating CK
                CK(i, i) += U(ir, 0) * U(ir, 0);
            }
        }

        CK(i+L,i+L)=1-CK(i,i);
    }





    //Calculating correlators C_{m,i} C_{i,m} i not Majorana,
    for (int i=0; i<L;++i){

        CK(i,2*L)=0;
        //  a^\dagger c_1   When there is NO Majorana
        for(int ir=0; ir<Lstate; ++ir){
            Decima_to_binary(Vect(ir), State,L+1);
            int count = 1;
            if (State(L)==0 &&  State(i)==1){

                State(i)=0;
                //Checking sign
                for(int ii=i; ii<L; ++ii){
                    if (State(ii)==1) {
                        count *= -1;
                    }
                }
                State(L)=1;

                double Vector=Binary_to_decimal(State, L+1);

                // Calculating CK
                for(int il=0; il<Lstate;++il){
                    if (Vector==Vect(il)){

                        CK(i,2*L)+=count*U(il,0)*U(ir,0);
                        break;
                    }
                }
            }
        }

        CK(i,2*L+1)=0;
        //  a c_1    When there is Majorana
        for(int ir=0; ir<Lstate; ++ir){
            Decima_to_binary(Vect(ir), State,L+1);
            int count = 1;
            if (State(L)==1 &&  State(i)==1){

                State(i)=0;
                //Checking sign
                for(int ii=i; ii<L+1; ++ii){
                    if (State(ii)==1) {
                        count *= -1;
                    }
                }
                State(L)=0;

                double Vector=Binary_to_decimal(State, L+1);

                // Calculating CK
                for(int il=0; il<Lstate;++il){
                    if (Vector==Vect(il)){

                        CK(i,2*L+1)+=count*U(il,0)*U(ir,0);
                        break;
                    }
                }
            }
        }

        CK(2*L,i)=CK(i,2*L);
        CK(2*L+1,i)=CK(i,2*L+1);

        CK(i+L,2*L)=-CK(i,2*L);
        CK(i+L,2*L+1)=-CK(i,2*L+1);
        CK(2*L,i+L)=CK(i+L,2*L);
        CK(2*L+1,i+L)=CK(i+L,2*L+1);

    }


    //Calculating correlators K(i,j) i,j not Majorana
    for (int i=0; i<L;++i){
        for(int j=i; j<L;++j){

            CK(i,j+L)=0;

            for(int ir=0; ir<Lstate; ++ir){
                Decima_to_binary(Vect(ir), State,L+1);
                int count = 1;
                if (State(j)==1 &&  State(i)==1){

                    //Checking sign
                    for(int ii=i; ii<j; ++ii){
                        if (State(ii)==1) {
                            count *= -1;
                        }
                    }
                    State(j)=0;
                    State(i)=0;
                    double Vector=Binary_to_decimal(State, L+1);

                    // Calculating CK
                    for(int il=0; il<Lstate;++il){
                        if (Vector==Vect(il)){

                            CK(i,j+L)+=count*U(il,0)*U(ir,0);
                            break;
                        }
                    }
                }
            }

            CK(j,i+L)=-CK(i,j+L);

            CK(j+L,i)=CK(i,j+L);
            CK(i+L,j)=CK(j,i+L);
        }
    }



}

void Decima_to_binary(double Vect, Eigen::VectorXi& State, int LL){
        int n=Vect;
        int count=0;
        for(int i=0; i<LL; i++){
            State(i) = n%2 ;
            n= n/2;
        }
}

double Binary_to_decimal(Eigen::VectorXi& State, int LL){
    int Vector=0;
    for (int i=0; i<LL;++i){
        Vector +=State(i)* pow(2,i);
    }

    return Vector;
}
