#include <iostream>
#include <fstream>
#include <unistd.h>
#include <fcntl.h>
#include <cmath>
#include <iomanip>
#include <sys/time.h>
using namespace std;

double f(double t, double x)
{
return 0.0;
//return sin(M_PI * x) * cos(3 * M_PI * x) * 2 * t + (1 + t * t) * (10 * M_PI * M_PI * sin(M_PI * x) * cos(3 * M_PI * x) + 6 * M_PI * M_PI * cos(M_PI * x) * sin(3 * M_PI * x));
//return (sin(t) * (x * x * x - x * x + 1) - (2 - cos(t)) * (6 * x - 2));
}

double u_0(double x)
{
//return exp(x);
//return sin(M_PI * x) * cos(3 * M_PI * x);
return x * x * x - x * x;
}

double phi_1(double t)
{
return 0.0;
//return 2 - cos(t);
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
return M_PI * (1 + t * t);
//return 2 - cos(t);
}

class Matrix {
    int n;
    double ** values;

public:
    double * optim;
	double ** X;

    Matrix(const int m = 0) {
        n = m;
        values = new double * [n];
		X = new double * [n];
        optim = new double [n];

        for (int i = 0; i < n; ++i) {

            values[i] = new double [n];
			X[i] = new double [n];

			for (int j = 0; j < n; ++j) {
				X[i][j] = 0;
			}
			X[i][i] = 1;
        }
    }

    void set_values();
	void set_optim();
	void set_ij(int * i, int * j, int mode);
	void lambda_appr(double sin, double cos, int i, int j);
	void eigenvectors(double sin, double cos, int i, int j);
	void eigenvalues(double eps, int mode);
    int find_max();

    ~Matrix() {
        for (int i = 0; i < n; ++i) {
            delete[] values[i];
			delete[] X[i];
        }
		delete[] X;
        delete[] optim;
        delete[] values;
    }
};

void Matrix::set_values() {
	for (int i = 0; i < n; ++i) {
		values[0][i] = 0;
		values[n - 1][i] = 0;
	}
	values[0][0] = -2e-4;
	values[0][1] = 1e-4;
	values[n - 1][n - 1] = -2e-4;
	values[n - 1][n - 2] = 1e-4;
	for (int i = 1; i < n - 1; ++i) {
		for (int j = 0; j < n; ++j) {
			values[i][j] = 0;
		}
		values[i][i - 1] = 1e-4;
		values[i][i] = -2e-4;
		values[i][i + 1] = 1e-4;
	}
}

void Matrix::set_optim() {
	for (int i = 0; i < n; ++i) {
		register double temp = 0;
		for (int j = 0; j < n; ++j) {
			if (j == i) continue;
			temp += values[i][j]*values[i][j];
		}
		optim[i] = temp;
	}
}

void Matrix::lambda_appr(double sin, double cos, int i, int j) {
	double temp = 0.0;
		temp = values[i][i];
		values[i][i] = values[i][i]*cos - values[j][i]*sin;
		values[j][i] = values[j][i]*cos + temp*sin;
		temp = values[i][j];
		values[i][j] = values[i][j]*cos - values[j][j]*sin;
		values[j][j] = values[j][j]*cos + temp*sin;

		temp = values[i][i];
		values[i][i] = values[i][i]*cos - values[i][j]*sin;
		values[i][j] = values[i][j]*cos + temp*sin;
		temp = values[j][i];
		values[j][i] = values[j][i]*cos - values[j][j]*sin;
		values[j][j] = values[j][j]*cos + temp*sin;

		for (int k = 0; k < n; ++k) {
			if ((k == i) || (k == j)) continue;

			temp = values[k][i];
			values[k][i] = values[k][i]*cos - values[k][j]*sin;
			values[k][j] = temp*sin + values[k][j]*cos;
			temp = values[i][k];
			values[i][k] = values[i][k]*cos - values[j][k]*sin;
			values[j][k] = temp*sin + values[j][k]*cos;

		}

	for (int k = 0; k < n; ++k) {
			if (k == i) {
				optim[k] += values[k][j]*values[k][j];
				continue;
			}

			if (k == j) {
				optim[k] += values[k][i]*values[k][i];
				continue;
			}

			optim[k] += (values[k][i]*values[k][i] + values[k][j]*values[k][j]);
		}

	for (int k = 0; k < n; ++k) {
			if (k == i) {
				optim[j] += values[j][k]*values[j][k];
				continue;
			}

			if (k == j) {
				optim[i] += values[i][k]*values[i][k];
				continue;
			}

			optim[j] += values[j][k]*values[j][k];
			optim[i] += values[i][k]*values[i][k];
		}

}

void Matrix::eigenvectors(double sin, double cos, int i, int j) {
	double temp = 0.0;
	for (int k = 0; k < n; ++k) {
		temp = X[k][i];
		X[k][i] = X[k][i]*cos - X[k][j]*sin;
		X[k][j] = temp*sin + X[k][j]*cos;
	}
}

void Matrix::set_ij(int * i_ptr, int * j_ptr, int mode) {

	int i = *i_ptr, j = *j_ptr;

	if (mode == 1) {
		double max = abs(values[i][j]);
		for (int k = 0; k < n; ++k) {
			for (int l = 0; l < n; ++l) {
				if (k == l) continue;
				if (abs(values[k][l]) > max) {
					max = abs(values[k][l]);
					i = k;
					j = l;
				}
			}
		}
	}

	if (mode == 2) {
		do {
			if (++j == n) {
				if (++i == n - 1) {
					i = 0;
					j = 1;
				}
				else j = i + 1;
			}
		}
		while (values[i][j] == 0.0);
	}

	if (mode == 3) {
		double max = optim[0];
		int max_ind = 0;

		for (int k = 1; k < n; ++k) {
			if (max < optim[k]) {
				max = optim[k];
				max_ind = k;
			}
		}
		i = max_ind;

		max = values[i][0];
		max_ind = 0;
		for (int k = 1; k < n; ++k) {
			if (k == i) continue;

			if (max < abs(values[i][k])) {
				max = abs(values[i][k]);
				max_ind = k;
			}
		}
		j = max_ind;
	}

	*i_ptr = i;
	*j_ptr = j;
}

void Matrix::eigenvalues(double eps, int mode) {
	int i = 0, j = 1, iter_cnt = 0;
	double x = 0.0, y = 0.0, accuracy = 0.0, sin = 0.0, cos = 0.0;

	struct timeval start, stop;
	gettimeofday(&start, NULL);

	for (int k = 0; k < n; ++k) {
		accuracy += optim[k];
	}

	while (accuracy >= eps) {
		set_ij(&i, &j, mode);

			for (int k = 0; k < n; ++k) {
				if (k == i) {
					optim[k] -= values[k][j]*values[k][j];
					continue;
				}

				if (k == j) {
					optim[k] -= values[k][i]*values[k][i];
					continue;
				}

				optim[k] -= (values[k][i]*values[k][i] + values[k][j]*values[k][j]);
			}

			optim[i] = 0;
			optim[j] = 0;


		x = - 2*values[i][j];
		y = values[i][i] - values[j][j];

		if (y == 0) {
			cos = (double)sqrt(2) / 2;
			sin = (double)sqrt(2) / 2;
		}
		else {
			register double temp = (double)2 * sqrt (x*x + y*y);

			cos = (double)sqrt((0.5) + (double)abs(y) / temp);

			if (((x >= 0) && (y >= 0)) || ((x < 0) && (y < 0)))
					   sin = abs(x) / (cos * temp);
			else sin = - abs(x) / (cos * temp);
		}

		lambda_appr(sin, cos, i, j);

		eigenvectors(sin, cos, i, j);

		accuracy = 0.0;
		for (int k = 0; k < n; ++k) {

			accuracy += optim[k];
		}

		++iter_cnt;
	}
	gettimeofday(&stop, NULL);

	cout << "Number of iterations on eigenvalues:" << iter_cnt << endl << endl;

	cout << "Execution time on eigenvalues:" << (((stop.tv_sec - start.tv_sec)*1000000) + stop.tv_usec - start.tv_usec)/(1.0*1000000) << "s" << endl << endl;
}

int Matrix::find_max() {
	double min = abs(values[0][0]);
	int min_ind = 0;
	for (int i = 0; i < n; ++i) {
		if (min > abs(values[i][i])) {
			min = abs(values[i][i]);
			min_ind = i;
		}
	}

	return min_ind;
}

void visualization_vect(double *vect, int type_left, int type_right, int N, int M, int mode, string S) {
	double *v_prev = new double[M + 1];
	double *v = new double[M + 1];
	double tau = 1/(N*1.0), h = 1/(M*1.0), hi1 = 0.0, mu1 = 0.0, hi2 = 0.0, mu2 = 0.0;

	double *alpha = new double[M + 1], *beta = new double[M + 1], Ci = (h*h + 2*tau), Ai = tau, Bi = tau;
	double *F = new double[M];

	ofstream fout;
	fout.open(S);

	for (int m = 1; m < M; ++m) {
		v_prev[m] = vect[m];
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

		register double norm_v = 0.0;
		for (int m = 0; m <= M; ++m){
			v_prev[m] = v[m];
			norm_v += v[m]*v[m];
		}
		fout << n*tau << "  " << sqrt(norm_v) << endl;
	}

	delete[] v;
	delete[] v_prev;
	delete[] alpha;
	delete[] beta;
	delete[] F;
	fout.close();
}

double *erase_domin(double *v_0, double *ksi_1, int n) {

	double numer = 0.0, denom = 0.0;
	for (int i = 0; i < n; ++i) {
		numer += v_0[i] * ksi_1[i];
		denom += ksi_1[i] * ksi_1[i];
	}
	double C1 = numer/(1.0*denom);

	double *w_0 = new double[n];
	for (int i = 0; i < n; ++i) {
		w_0[i] = v_0[i] - C1*ksi_1[i];
	}
	return w_0;
}

int main(int argc, char **argv) {

/*Command line format: ./bin_file key N M*/

	if (argc < 2) {
		perror("No enough arguments");
		return -1;
	}

	int N = atoi(argv[1]), M = atoi(argv[2]);

	double *v = new double[M + 1], *ksi_1 = new double[M + 1], h = 1/(1.0*M);

	for (int i = 0; i <= M; ++i) {
		v[i] = u_0(i*h);
	}
	visualization_vect(v, 1, 1, N, M, 2, "data.txt");

	Matrix A(M - 1);
	A.set_values();
	A.set_optim();
	A.eigenvalues(1e-24, 3);

	int ind = A.find_max();
	ksi_1[0] = 0;
	ksi_1[M] = 0;
	for (int i = 1; i < M; ++i) {
		ksi_1[i] = A.X[i - 1][ind];
	}
	double *w = erase_domin(v, ksi_1, M + 1);
	visualization_vect(w, 1, 1, N, M, 2, "datav.txt");

	delete[] v;
	delete[] ksi_1;
	delete[] w;
}
