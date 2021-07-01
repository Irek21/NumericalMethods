#include <iostream>
#include <fstream>
#include <unistd.h>
#include <csignal>
#include <cmath>
using namespace std;

double u(double t, double x)
{
return (1 + (t * t)) * sin(M_PI * x) * cos(3 * M_PI * x);
//return (2 - cos(t)) * (x * x * x - x * x + 1);
}

double f(double t, double x)
{
//return sin(M_PI * x) * cos(3 * M_PI * x) * 2 * t + (1 + t * t) * (10 * M_PI * M_PI * sin(M_PI * x) * cos(3 * M_PI * x) + 6 * M_PI * M_PI * cos(M_PI * x) * sin(3 * M_PI * x));
//return (sin(t) * (x * x * x - x * x + 1) - (2 - cos(t)) * (6 * x - 2));
	return 0;
}

double u_0(double x)
{
//return sin(M_PI * x) * cos(3 * M_PI * x);
//return x * x * x - x * x + 1;
	return 0;
}

double phi_1(double t)
{
//return 0.0;
//return 2 - cos(t);
	return sin(t);
}

double phi_2(double t)
{
return 0.0;
//return 2 - cos(t);
}

double psi_1(double t)
{
return M_PI * (1 + t * t);
//return 0.0;
}

double psi_2(double t)
{
//return M_PI * (1 + t * t);
//return 2 - cos(t);
	return 0;
}

void visualization(int type_left, int type_right, int N, int M, int mode) {
	double *v_prev = new double[M + 1];
	double *v = new double[M + 1];
	double tau = 1/(N*1.0), h = 1/(M*1.0), hi1 = 0.0, mu1 = 0.0, hi2 = 0.0, mu2 = 0.0;
	
	double *alpha = new double[M + 1], *beta = new double[M + 1], Ci = (h*h + 2*tau), Ai = tau, Bi = tau;
	double *F = new double[M];
	
	ofstream fout;
	fout.open("data.txt");
	
	for (int m = 1; m < M; ++m) {
		v_prev[m] = u_0(m*h);
	}
	
	if (mode == 1) {
		v_prev[0] = (type_left == 1) ? phi_1(0) : (2*psi_1(0)*h + v_prev[2] - 4*v_prev[1])/(-3.0);
		v_prev[M] = (type_right == 1) ? phi_2(0) : (2*psi_2(0)*h - v_prev[M - 2] + 4*v_prev[M - 1])/(3.0);
	} else {
		v_prev[0] = (type_left == 1) ? phi_1(0) : v_prev[1] - h*psi_1(0);
		v_prev[M] = (type_right == 1) ? phi_2(0) : v_prev[M - 1] - h*psi_2(0);
	}
	
	for (int n = 1; n <= 10*N; ++n){
		
		switch (mode) {
			case 1: {
				for (int m = 1; m < M; ++m) {
					v[m] = v_prev[m] + tau*((v_prev[m + 1] - 2*v_prev[m] + v_prev[m - 1])/(h*h*1.0) + f(n*tau, m*h));
				}
				v[0] = (type_left == 1) ? phi_1(n*tau) : (2*psi_1(n*tau)*h + v[2] - 4*v[1])/(-3.0);
				v[M] = (type_right == 1) ? phi_2(n*tau) : (2*psi_2(n*tau)*h - v[M - 2] + 4*v[M - 1])/(3.0);
				
				break;
			}
			
			case 2: {
				for (int m = 1; m <= M - 1; ++m) {
					F[m] = h*h*(f(n*tau, m*h)*tau + v_prev[m]);
				}
				
				if (type_left == 1) {
					hi1 = 0;
					mu1 = phi_1(n*tau);
				} else {
					hi1 = 1;
					mu1 = -h*psi_1(n*tau);
				}
				
				if (type_right == 1) {
					hi2 = 0;
					mu2 = phi_2(n*tau);
				} else {
					hi2 = 1;
					mu2 = h*psi_2(n*tau);
				}
			}
				
			alpha[1] = hi1;		
			beta[1] = mu1;

			for (int m = 1; m <= M - 1; ++m) {		
				register double denom = - Ai*alpha[m] + Ci;			
				alpha[m + 1] = Bi/(1.0*denom);
				beta[m + 1] = (F[m] + Ai*beta[m])/(1.0*denom);

			}

			v[M] = (mu2 + hi2*beta[M])/(1.0*(1 - hi2*alpha[M]));
			for (int m = M - 1; m >= 0; --m) {			
				v[m] = v[m + 1]*alpha[m + 1] + beta[m + 1];		
			}
		}
		
		for (int m = 0; m <= M; ++m){
			v_prev[m] = v[m];
			fout << m*h << "  " << v[m] << endl; 
		}
		
		/*if (n == 50) {
			char c;
			fout.close();
			cout << "Program waits for a key..." << endl;
			cin >> c;
			fout.open("data.txt");
		}*/
		
		fout << endl << endl;
	}
	
	delete[] v;
	delete[] v_prev;
	delete[] alpha;
	delete[] beta;
	delete[] F;
	fout.close();
}

void testing(int type_left, int type_right, int N, int M, int mode) {
	double *v_prev = new double[M + 1];
	double *v = new double[M + 1];
	double tau = 1/(N*1.0), h = 1/(M*1.0), hi1 = 0.0, mu1 = 0.0, hi2 = 0.0, mu2 = 0.0;
	
	double *alpha = new double[M + 1], *beta = new double[M + 1], Ci = (h*h + 2*tau), Ai = tau, Bi = tau;
	double *F = new double[M];
	
	for (int m = 1; m < M; ++m) {
		v_prev[m] = u_0(m*h);
	}
	
	if (mode == 1) {
		v_prev[0] = (type_left == 1) ? phi_1(0) : (2*psi_1(0)*h + v_prev[2] - 4*v_prev[1])/(-3.0);
		v_prev[M] = (type_right == 1) ? phi_2(0) : (2*psi_2(0)*h - v_prev[M - 2] + 4*v_prev[M - 1])/(3.0);
	} else {
		v_prev[0] = (type_left == 1) ? phi_1(0) : v_prev[1] - h*psi_1(0);
		v_prev[M] = (type_right == 1) ? phi_2(0) : v_prev[M - 1] - h*psi_2(0);
	}
	
	for (int n = 1; n <= N; ++n){
		
		switch (mode) {
			case 1: {
				for (int m = 1; m < M; ++m) {
					v[m] = v_prev[m] + tau*((v_prev[m + 1] - 2*v_prev[m] + v_prev[m - 1])/(h*h*1.0) + f(n*tau, m*h));
				}
				v[0] = (type_left == 1) ? phi_1(n*tau) : (2*psi_1(n*tau)*h + v[2] - 4*v[1])/(-3.0);
				v[M] = (type_right == 1) ? phi_2(n*tau) : (2*psi_2(n*tau)*h - v[M - 2] + 4*v[M - 1])/(3.0);
				
				break;
			}
			
			case 2: {
				for (int m = 1; m <= M - 1; ++m) {
					F[m] = h*h*(f(n*tau, m*h)*tau + v_prev[m]);
				}
				
				if (type_left == 1) {
					hi1 = 0;
					mu1 = phi_1(n*tau);
				} else {
					hi1 = 1;
					mu1 = -h*psi_1(n*tau);
				}
				
				if (type_right == 1) {
					hi2 = 0;
					mu2 = phi_2(n*tau);
				} else {
					hi2 = 1;
					mu2 = h*psi_2(n*tau);
				}
			}
				
			alpha[1] = hi1;		
			beta[1] = mu1;

			for (int m = 1; m <= M - 1; ++m) {		
				register double denom = - Ai*alpha[m] + Ci;			
				alpha[m + 1] = Bi/(1.0*denom);
				beta[m + 1] = (F[m] + Ai*beta[m])/(1.0*denom);

			}

			v[M] = (mu2 + hi2*beta[M])/(1.0*(1 - hi2*alpha[M]));
			for (int m = M - 1; m >= 0; --m) {			
				v[m] = v[m + 1]*alpha[m + 1] + beta[m + 1];		
			}
		}
		
		if (n == N) continue;
		for (int m = 0; m <= M; ++m){
			v_prev[m] = v[m];
		}
	}

	double sub = 0.0, C_h = 0.0, sum_sq = 0.0, l_2h = 0.0;
	double sum_reg = 0.0, C_h_reg = 0.0, l_2h_reg;
	for (int m = 0; m <= M; ++m){
		register double func = u(1, m*h);
		if ((sub = abs(func - v[m])) > C_h) C_h = sub;
		sum_sq += sub*sub;
		if (abs(func) > C_h_reg) C_h_reg = abs(func);
		sum_reg += func*func;
	}
	l_2h = sqrt(h*sum_sq*1.0);
	l_2h_reg = sqrt(h*sum_reg*1.0);
	
	cout << "Ch absolute: " << C_h << endl;
	cout << "l2h absolute: " << l_2h << endl;
	cout << "Ch regarding: " << C_h/(1.0*C_h_reg) << endl;
	cout << "l2h regarding: " << l_2h/(1.0*l_2h_reg) << endl;
	//cout << v[M/2 + 1] << endl;
	
	delete[] v;
	delete[] v_prev;
	delete[] alpha;
	delete[] beta;
	delete[] F;
}

int main(int argc, char ** argv) 
{
	if (argc < 7) {
		perror("No enough arguements");
		return -1;
	}
	int type_left = atoi(argv[1]), type_right = atoi(argv[2]), N = atoi(argv[3]), M = atoi(argv[4]), mode = atoi(argv[5]), regime = atoi(argv[6]);
	if (regime == 1) {
		testing(type_left, type_right, N, M, mode);
	} else visualization(type_left, type_right, N, M, mode);
}