
// g++ mpc_test.cpp -o mpctest -std=c++11 -larmadillo

#include <armadillo>
#include <iostream>
#include <stdlib.h>

using namespace std;
using namespace arma;

mat A_d_all;
mat B_d;

mat H, C;

int HORIZON;
float dT = 0.01;

void blkdiag(mat, mat, mat&);
void bldMatrixDiag(float*, mat&, int);
void pushMatrix_row(mat&, mat);
void pushMatrix_col(mat&, mat);
// void delete_arr(float*);

int main () {

        // settings
        dT = 0.01;
        HORIZON = 10;

        //assign matrices
        mat A_d_triple = {{1, dT, 0}, {0,1,0}, {dT, dT*dT/2, 1}};

        //assign limitations
        float u1_min = -10;
        float u1_max = 10;
        float u2_min = -10;
        float u2_max = 10;
        float u3_min = -5;
        float u3_max = 20;

        float z_max = -15;
        float z_min = -25;

        float z_d = -20;
        float mu = 0.01;

        float tan_half_phi = 0.5;
        float tan_half_psi = 0.5;

        mat tmp;
        blkdiag( A_d_triple, A_d_triple, tmp );
        blkdiag( A_d_triple, tmp, A_d_all);

        // A_d_all.print("a-d-all");

        mat B_one_state;
        B_one_state << dT*dT/2 <<endr
                                << dT <<endr
                                << dT*dT*dT/6 <<endr;

        blkdiag( B_one_state, B_one_state, tmp );
        blkdiag( B_one_state, tmp, B_d);

        mat R(3,3);
        R.eye();
        R = R*0.001;

        float a[] = {1,1,0.001,1,1,0.001,0.001,0.01,0.001};

        mat Q_end, Q;
        bldMatrixDiag(a, Q_end, sizeof(a)/sizeof(*a));

        float b[] = {1, 0.05, 0.01, 1, 0.05, 0.01, 0.01, 0.01, 0.01};
        // cout <<"size" << sizeof(b)/sizeof(*b)<<endl;

        bldMatrixDiag(b, Q, sizeof(b)/sizeof(*b));

        H.reset();
        C.reset();
        mat zero_mat, eye_mat;

        for (int i = 0; i < HORIZON; i ++) {

                if(i != HORIZON - 1) {

                        blkdiag(H, R, H);
                        blkdiag(H, Q, H);
                        H.print("H");

                } else {

                        blkdiag(H, R, H);
                        blkdiag(H, Q_end, H);

                }

                if ( i == 0) {

                        tmp.reset();
                        pushMatrix_row(tmp, -B_d);

                        eye_mat.eye(9,9);
                        eye_mat.print("eye_mat");

                        cout<<tmp.n_rows<<endl;

                        pushMatrix_row(tmp, eye_mat);

                        tmp.print("aa");
                        cout<<__LINE__<<endl;

                        eye_mat.zeros(9, 12*(HORIZON - 2));

                        pushMatrix_row(tmp, eye_mat);
                        cout<<__LINE__<<endl;

                        pushMatrix_col(C, tmp);
                        cout<<__LINE__<<endl;


                }

                else if (i == HORIZON - 1) {

                        tmp.reset();
                        zero_mat.zeros(9,3);

                        pushMatrix_row(tmp, zero_mat);
                        zero_mat.zeros(9, 12*(i - 2));
                        pushMatrix_row(tmp, zero_mat);
                        pushMatrix_row(tmp, -A_d_all);
                        pushMatrix_row(tmp, -B_d);

                        eye_mat.eye(9,9);
                        pushMatrix_row(tmp, eye_mat);

                        pushMatrix_col(C, tmp);

                }

                else if (i == 1) {

                        cout<<__LINE__<<endl;

                        tmp.reset();
                        zero_mat.zeros(9, 3);
                        cout<<__LINE__<<endl;

                        pushMatrix_row(tmp, zero_mat);
                        cout<<__LINE__<<endl;

                        pushMatrix_row(tmp, -A_d_all);
                        pushMatrix_row(tmp, -B_d);

                        eye_mat.eye(9,9);
                        pushMatrix_row(tmp, eye_mat);
                        zero_mat.zeros(9, 12 * (HORIZON - 2));
                        pushMatrix_row(tmp, zero_mat);

                        pushMatrix_col(C, tmp);

                }

                else {

                        tmp.reset();
                        zero_mat.zeros(9, 3);
                        pushMatrix_row(tmp, zero_mat);
                        zero_mat.zeros(9, (i - 2)*12);
                        pushMatrix_row(tmp, zero_mat);
                        pushMatrix_row(tmp, -A_d_all);
                        pushMatrix_row(tmp, -B_d);
                        eye_mat.eye(9, 9);
                        pushMatrix_row(tmp, eye_mat);
                        zero_mat.zeros(9, 12*(HORIZON - i));
                        pushMatrix_row(tmp, zero_mat);

                        pushMatrix_col(C, tmp);

                }


        }

}

void blkdiag(mat A, mat B, mat& C) {

        mat C_tmp;
        C_tmp.resize(A.n_rows + B.n_rows, A.n_cols + B.n_cols);
        C_tmp.zeros();

        for( int i = 0; i < A.n_rows; i ++) {

                for (int j = 0; j < A.n_cols; j ++) {

                        C_tmp.at(i,j) = A.at(i,j);

                }

        }

        for( int i = 0; i < B.n_rows; i ++) {

                for( int j = 0; j < B.n_cols; j ++ ) {

                        C_tmp.at( i + A.n_rows, j + A.n_cols ) = B.at(i,j);

                }

        }

        C = C_tmp;

}

void bldMatrixDiag(float *Vec, mat& MATRIX, int Length) {

        MATRIX.zeros(Length, Length);

        for (int i = 0; i < Length; i ++) {

                MATRIX.at(i,i) = Vec[i];

        }

}

void pushMatrix_row(mat& MATRIX, mat ToBePushed) {


        cout<< "MATRIX" << MATRIX.n_rows<< ", " <<MATRIX.n_cols <<endl;
        cout<< "ToBePushed" << ToBePushed.n_rows << ", " << ToBePushed.n_cols <<endl;

        if (MATRIX.is_empty()) {

                cout <<__LINE__<<endl;
                MATRIX.resize(ToBePushed.n_rows, ToBePushed.n_cols);
                for(int i = 0; i < ToBePushed.n_rows; i ++) {

                        for (int j = 0; j < ToBePushed.n_cols; j ++) {

                                MATRIX.at(i, j) = ToBePushed.at(i, j);

                        }

                }

                return;

        }

        if (MATRIX.n_rows != ToBePushed.n_rows) {

                cout << " wrong dimension!" <<endl;
                abort();

        }

        mat tmp;
        tmp.zeros(MATRIX.n_rows, MATRIX.n_cols + ToBePushed.n_cols);

        cout<<"tmp: " << tmp.n_rows <<", " <<tmp.n_cols <<endl;
        for(int i = 0; i < MATRIX.n_rows; i ++) {

                cout << i <<endl;
                for(int j = 0; j < MATRIX.n_cols; j ++) {

                        tmp.at(i,j) = MATRIX.at(i,j);

                }

        }

        for(int i = 0; i < ToBePushed.n_rows; i  ++) {

                for(int j = 0; j < ToBePushed.n_cols; i ++) {

                        tmp.at(i, j + MATRIX.n_cols) = ToBePushed.at(i, j);
                                                cout<<__LINE__<<endl;

                }

        }

        MATRIX = tmp;

}

void pushMatrix_col(mat& MATRIX, mat ToBePushed) {


        if (MATRIX.is_empty()) {

                MATRIX.resize(ToBePushed.n_rows, ToBePushed.n_cols);
                for(int i = 0; i < ToBePushed.n_rows; i ++) {

                        for (int j = 0; j < ToBePushed.n_cols; j ++) {

                                MATRIX.at(i, j) = ToBePushed.at(i, j);

                        }

                }

                return;

        }

        if(MATRIX.n_cols != ToBePushed.n_cols) {

                cout << " wrong dimension!" <<endl;
                abort();

        }

        mat tmp;

        for (int i = 0; i < MATRIX.n_rows; i ++) {

                for (int j = 0; j < MATRIX.n_cols; j ++) {

                        tmp.at(i,j) = MATRIX.at(i,j);

                }

        }

        for (int i = 0; i < ToBePushed.n_rows; i ++) {

                for ( int j = 0; j < ToBePushed.n_rows; j ++) {

                        tmp.at(i, j + MATRIX.n_cols) = ToBePushed.at(i,j);

                }

        }

        MATRIX = tmp;


}
// void delete_arr(float* Vec) {

// 	int Length = sizeof(Vec)/sizeof(*Vec);
// 	for(int i = 0; i < Length; i ++) {

// 		free(Vec[i]);

// 	}

// }
