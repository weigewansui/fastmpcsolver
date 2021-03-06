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
void initialize_Matrices();
// void delete_arr(float*);
// 

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

// Coefficients for the barrier

float mu_state = 0.1;
float mu_u = 0.0001;

int main() {

	initialize_Matrices();
	cout << C.n_rows<< ", " <<C.n_cols<<endl;

	// initial conditions
	mat tmp0;
	tmp0 = {0, 0, 0, 1, 0, 0, 1, 0, 0, -16, 0, 0};
	tmp0 = trans(tmp0);

	mat z0(12*HORIZON, 1);
	pushMatrix_col(z0, tmp0);

	mat beq(9*HORIZON,1);

	mat x0 = {1, 0, 0, 1, 0, 0, -16, 0, 0};
	x0 = trans(x0);

	//for testing purpose
	float noise_mu = 0.0001;
	mat noise = noise_mu * x0;
	mat beqFirstHorizon = A_d_all*x0 + noise;
	
	for(int i = 0; i < 9; i++) {

		beq(i,0) = beqFirstHorizon(i, 0);
		z0(i,0) = tmp0(i,0);

	}

	// set the values of entrices of beq for other horizons
	for(int i = 9; i < HORIZON*9; i ++) {

		beq(i, 0) = noise(i%9,0);
		z0(i,0) = tmp0(i%9,0);

	}

	mat v0 =  0.05 * randi<mat>(9*HORIZON, 1, distr_param(0, 20));

	mat z = z0;
	mat v = v0;

 	vec r(21);
 	r.ones();
	vec r_d(12*HORIZON);
	vec r_p(9*HORIZON);

	float w_r_p = 1e-4;
	float w_d_p = 1e-4;

	int loop_num = 0;

	float fval;
	float x, y, z_h;
	mat HESSIAN, GRADIENT;

	float u1, u2, u3, h;
	float tmp_denum_grad_x_1, tmp_denum_grad_x_2, tmp_denum_grad_y_1, tmp_denum_grad_y_2;
	float d_phi_dx, d_phi_dy, d_phi_dz;
	float d_phi_du1, d_phi_du2, d_phi_du3;

	float tmp_denum_x_1, tmp_denum_x_2, tmp_denum_y_1, tmp_denum_y_2;
	float dd_phi_dxdx, dd_phi_dxdz, dd_phi_dydy, dd_phi_dydz, dd_phi_dzdz;
	float dd_phi_du1du1, dd_phi_du2du2, dd_phi_du3du3;

	while(norm(r) > 1e-4) {

		// cout <<__LINE__<<endl;
		fval = dot(trans(z), (H * z));
		// cout <<__LINE__<<endl;


		HESSIAN.zeros(12*HORIZON, 12*HORIZON);
		GRADIENT.zeros(12*HORIZON, 1);
		// cout <<__LINE__<<endl;


		for (int i = 0; i < HORIZON; i++) {

			x = z.at(3 + i * 12);
		// cout <<__LINE__<<endl;

			y = z.at(6 + i * 12);
		// cout <<__LINE__<<endl;

			z_h = z (9 + i * 12);
		// cout <<__LINE__<<endl;

			u1 = z.at(0 + 12 * i);
		// cout <<__LINE__<<endl;


			u2 = z.at(1 + 12 * i);
		// cout <<__LINE__<<endl;

			u3 = z.at(2 + 12 * i);
		// cout <<__LINE__<<endl;


			h = -z_h;

			tmp_denum_grad_x_1 = h*tan_half_phi - x;
			tmp_denum_grad_x_2 = h*tan_half_phi + x;
			tmp_denum_grad_y_1 = h*tan_half_psi - y;
			tmp_denum_grad_y_2 = h*tan_half_psi + y;

			d_phi_dx = 1/tmp_denum_grad_x_1 - 1/tmp_denum_grad_x_2;
         	d_phi_dy = 1/tmp_denum_grad_y_1 - 1/tmp_denum_grad_y_2;
         	d_phi_dz = 1/(z_max - z_h) - 1/(z_h - z_min) + tan_half_phi/tmp_denum_grad_x_1 + tan_half_phi/tmp_denum_grad_x_2
             	+ tan_half_psi/tmp_denum_grad_y_1 + tan_half_psi/tmp_denum_grad_y_2;

         	d_phi_du1 = 1/(u1_max - u1) - 1/(u1 - u1_min);
        	d_phi_du2 = 1/(u2_max - u2) - 1/(u2 - u2_min);
         	d_phi_du3 = 1/(u3_max - u3) - 1/(u3 - u3_min);

         	GRADIENT.at(i*12) = mu_u*d_phi_du1;
         	GRADIENT.at(i*12 + 1) = mu_u*d_phi_du2;
         	GRADIENT.at(i*12 + 2) = mu_u*d_phi_du3;
         	GRADIENT.at(i*12 + 3) = mu_state*d_phi_dx;
         	GRADIENT.at(i*12 + 6) = mu_state*d_phi_dy;
         	GRADIENT.at(i*12 + 9) = mu_state*d_phi_dz;
		// cout <<__LINE__<<endl;


         	tmp_denum_x_1 = (h*tan_half_phi - x)*(h*tan_half_phi - x);
         	tmp_denum_x_2 = (h*tan_half_phi + x)*(h*tan_half_phi + x);

         	tmp_denum_y_1 = (h*tan_half_psi - y)*(h*tan_half_psi - y);
         	tmp_denum_y_2 = (h*tan_half_psi + y)*(h*tan_half_psi + y);

         	dd_phi_dxdx = 1/tmp_denum_x_1 + 1/tmp_denum_x_2;
         	dd_phi_dxdz = tan_half_phi/tmp_denum_x_1 - tan_half_phi/tmp_denum_x_2;

         	dd_phi_dydy = 1/tmp_denum_y_1 + 1/tmp_denum_y_2;
         	dd_phi_dydz = tan_half_psi/tmp_denum_y_1 - tan_half_psi/tmp_denum_y_2;

         	dd_phi_dzdz = 1/(z_max - z_h)/(z_max - z_h) + 1/(z_h - z_min)/(z_h - z_min) + tan_half_phi*tan_half_phi/tmp_denum_x_1
            	 + tan_half_phi*tan_half_phi/tmp_denum_x_2 + tan_half_psi*tan_half_psi/tmp_denum_y_1 + tan_half_psi*tan_half_psi/tmp_denum_y_2;

         	dd_phi_du1du1 = 1/(u1_max - u1)/(u1_max - u1) + 1/(u1 - u1_min)/(u1 - u1_min);
         	dd_phi_du2du2 = 1/(u2_max - u2)/(u2_max - u2) + 1/(u2 - u2_min)/(u2 - u2_min);
         	dd_phi_du3du3 = 1/(u3_max - u3)/(u3_max - u3) + 1/(u3 - u3_min)/(u3 - u3_min);

         	HESSIAN.at(12*i, 12*i) = mu_u*dd_phi_du1du1;
         	HESSIAN.at(12*i + 1, 12*i + 1) = mu_u*dd_phi_du2du2;
         	HESSIAN.at(12*i + 2, 12*i + 2) = mu_u*dd_phi_du3du3;

         	HESSIAN.at(12*i + 3, 12*i + 3) = mu_state*dd_phi_dxdx;
         	HESSIAN.at(12*i + 3, 12*i + 9) = mu_state*dd_phi_dxdz;
         	HESSIAN.at(12*i + 6, 12*i + 6) = mu_state*dd_phi_dydy;
         	HESSIAN.at(12*i + 6, 12*i + 9) = mu_state*dd_phi_dydz;
         	HESSIAN.at(12*i + 9, 12*i + 9) = mu_state*dd_phi_dzdz;
         	HESSIAN.at(12*i + 9, 12*i + 3) = mu_state*dd_phi_dxdz;
         	HESSIAN.at(12*i + 9, 12*i + 6) = mu_state*dd_phi_dydz;


         	// cout << HESSIAN.n_rows<< ", " <<HESSIAN.n_cols<<endl;
		}

		cout << z.n_rows<< ", " <<z.n_cols<<endl;
		cout << H.n_rows<< ", " <<H.n_cols<<endl;
		cout << C.n_rows<< ", " <<C.n_cols<<endl;
		cout << v.n_rows<< ", " <<v.n_cols<<endl;

		r_d = 2*H*z + GRADIENT + trans(C)*v;
		r_p = C*z - beq;

		r = r_d;
		r.insert_rows(12*HORIZON, r_p);

	}








}

void initialize_Matrices () {

	// settings
	dT = 0.01;
	HORIZON = 10;

	//assign matrices
	mat A_d_triple = {{1, dT, 0}, {0,1,0}, {dT, dT*dT/2, 1}};

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
			// H.print("H");

		} else {

			blkdiag(H, R, H);
			blkdiag(H, Q_end, H);

		}

		if ( i == 0) {

			tmp.reset();
			pushMatrix_row(tmp, -B_d);

			// tmp.print("tmp-B_d");
			
			eye_mat.eye(9,9);
			// eye_mat.print("eye_mat");

			// cout<<tmp.n_rows<<endl;

			pushMatrix_row(tmp, eye_mat);

			// tmp.print("aa");
			// cout<<__LINE__<<endl;

			eye_mat.zeros(9, 12*(HORIZON - 1));

			pushMatrix_row(tmp, eye_mat);
			// cout<<__LINE__<<endl;

			pushMatrix_col(C, tmp);
			// cout<<__LINE__<<endl;

			// cout << "tmp: " <<tmp.n_rows << ", " << tmp.n_cols <<endl;
			// cout << "C: " << C.n_rows << ", " <<C.n_cols << endl;

			// return 0;


		}

		else if (i == HORIZON - 1) {

			tmp.reset();
			zero_mat.zeros(9,3);

			pushMatrix_row(tmp, zero_mat);
			// cout << __LINE__ <<endl;

			zero_mat.zeros(9, 12*((i + 1) - 2));
			pushMatrix_row(tmp, zero_mat);
			pushMatrix_row(tmp, -A_d_all);
			pushMatrix_row(tmp, -B_d);

			// cout << __LINE__ <<endl;

			eye_mat.eye(9,9);
			pushMatrix_row(tmp, eye_mat);
			// cout << __LINE__ <<endl;

			// cout << "C: " <<C.n_rows <<", " <<C.n_cols <<endl;
			// cout << "tmp: " << tmp.n_rows <<", "<< tmp.n_cols<<endl;
			pushMatrix_col(C, tmp);	
			// cout << __LINE__ <<endl;


		}

		else if (i == 1) {

			// cout<<__LINE__<<endl;

			tmp.reset();
			zero_mat.zeros(9, 3);
			// cout<<__LINE__<<endl;

			pushMatrix_row(tmp, zero_mat);
			// cout<<__LINE__<<endl;

			// cout << "tmp: " <<tmp.n_rows << ", " << tmp.n_cols <<endl;
			// cout << "-A_d_all: " << A_d_all.n_rows << ", " <<A_d_all.n_cols << endl;
			pushMatrix_row(tmp, -A_d_all);
			pushMatrix_row(tmp, -B_d);

			// cout << __LINE__<<endl;

			eye_mat.eye(9,9);
			pushMatrix_row(tmp, eye_mat);
			// cout << __LINE__<<endl;
			
			zero_mat.zeros(9, 12 * (HORIZON - 2));
			pushMatrix_row(tmp, zero_mat);
			// cout << __LINE__<<endl;

			// cout << "tmp: " <<tmp.n_rows << ", " << tmp.n_cols <<endl;
			// cout << "C: " << C.n_rows << ", " <<C.n_cols << endl;
			// cout << __LINE__<<endl;

			pushMatrix_col(C, tmp);
			// cout << __LINE__<<endl;

		}

		else {

			tmp.reset();
			zero_mat.zeros(9, 3);
			pushMatrix_row(tmp, zero_mat);
			zero_mat.zeros(9, ((i + 1) - 2)*12);
			pushMatrix_row(tmp, zero_mat);
			// cout << __LINE__ <<endl;
			pushMatrix_row(tmp, -A_d_all);
			pushMatrix_row(tmp, -B_d);
			eye_mat.eye(9, 9);
			pushMatrix_row(tmp, eye_mat);
			zero_mat.zeros(9, 12*(HORIZON - (i + 1)));
			pushMatrix_row(tmp, zero_mat);

			pushMatrix_col(C, tmp);
			// cout << __LINE__ <<endl;

		}


	}

	cout << "C: " << C.n_rows << ", " <<C.n_cols << endl;

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


	// cout<< "MATRIX" << MATRIX.n_rows<< ", " <<MATRIX.n_cols <<endl;
	// cout<< "ToBePushed" << ToBePushed.n_rows << ", " << ToBePushed.n_cols <<endl;

	if (MATRIX.is_empty()) {

		// cout <<__LINE__<<endl;
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

	// cout<<"tmp: " << tmp.n_rows <<", " <<tmp.n_cols <<endl;
	for(int i = 0; i < MATRIX.n_rows; i ++) {

		// cout << i <<endl;
		for(int j = 0; j < MATRIX.n_cols; j ++) {

			tmp.at(i,j) = MATRIX.at(i,j);

		}

	}

	for(int i = 0; i < ToBePushed.n_rows; i  ++) {

		for(int j = 0; j < ToBePushed.n_cols; j ++) {

			tmp.at(i, j + MATRIX.n_cols) = ToBePushed.at(i, j);
						// cout<<__LINE__<<endl;

		}

	}

	MATRIX = tmp;

}

void pushMatrix_col(mat& MATRIX, mat ToBePushed) {

		// cout<<__LINE__<<endl;

	if (MATRIX.is_empty()) {

		MATRIX.resize(ToBePushed.n_rows, ToBePushed.n_cols);
		for(int i = 0; i < ToBePushed.n_rows; i ++) {

			for (int j = 0; j < ToBePushed.n_cols; j ++) {

				MATRIX.at(i, j) = ToBePushed.at(i, j);

			}

		}

		return;

	}

		// cout<<__LINE__<<endl;
	
	if(MATRIX.n_cols != ToBePushed.n_cols) {

		cout << " wrong dimension!" <<endl;
		abort();	

	}

	mat tmp;

	// cout<<__LINE__<<endl;
	tmp.zeros(MATRIX.n_rows + ToBePushed.n_rows, MATRIX.n_cols);
	
	// cout<<__LINE__<<endl;

	for (int i = 0; i < MATRIX.n_rows; i ++) {
		// cout<<__LINE__<<endl;

		for (int j = 0; j < MATRIX.n_cols; j ++) {

			tmp.at(i,j) = MATRIX.at(i,j);

		}

	}

	// cout<<__LINE__<<endl;
	// cout<<ToBePushed.n_rows<<", "<<ToBePushed.n_cols<<endl;
	// cout<<tmp.n_rows<<", " <<tmp.n_cols<<endl;

	for (int i = 0; i < ToBePushed.n_rows; i ++) {
	// cout<<__LINE__<<endl;


		for ( int j = 0; j < ToBePushed.n_rows; j ++) {

			tmp.at(i + MATRIX.n_rows, j) = ToBePushed.at(i,j);

		}

	}
	// cout<<__LINE__<<endl;
	// cout<<__LINE__<<endl;
	// cout<<ToBePushed.n_rows<<", "<<ToBePushed.n_cols<<endl;
	// cout<<tmp.n_rows<<", " <<tmp.n_cols<<endl;

	// MATRIX = tmp;
	// cout<<__LINE__<<endl;
	// cout<<__LINE__<<endl;


}