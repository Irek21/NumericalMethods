#include <iostream>
#include <unistd.h>
#include <fcntl.h>
#include <cmath>
#include <iomanip>
#include <sys/time.h>
using namespace std;

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

    void get_values();
    void set_values();
	void set_optim();
	void set_ij(int * i, int * j, int mode);
	void lambda_appr(double sin, double cos, int i, int j);
	void eigenvectors(double sin, double cos, int i, int j);
	void eigenvalues(double eps, int mode);
	void print_values();
    void print_eigen_sollutions();

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

double Function(int i, int j, int n) {
    return 1.0 / (2 * n - i - j -1);
}

void Matrix::set_values() {
	for (int i = 0; i < n - 1; ++i) {
		for (int j = 0; j < n; ++j) {
			values[i][j] = 0;
		}
		values[i][i] = 1;
		values[i][n - 1] = i + 1;
	}

	for (int i = 0; i < n; ++i) {
		values[n - 1][i] = i + 1;
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

void Matrix::get_values() {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cin >> values[i][j];
        }
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

	cout << "Number of iterations:" << iter_cnt << endl << endl;
	cout << "Execution time:" << (((stop.tv_sec - start.tv_sec)*1000000) + stop.tv_usec - start.tv_usec)/(1.0*1000000) << "s" << endl << endl;
}

void Matrix::print_eigen_sollutions() {
	cout << "Eigenvalues:" << endl;//setprecision(20) << endl;
	for (int i = 0; i < n; ++i) {
		cout << values[i][i] << "   ";
	}

	cout << endl;

	if (n > 10) {
		int fd_out = open("foutput.txt", O_RDWR | O_CREAT | O_TRUNC, 0777);
		dup2(fd_out, 1);
	}

	cout << "Eigenvectors:" << endl;
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			cout.width(10);
			cout << X[i][j] << "   ";
		}
		cout << endl << endl;
	}


}

void Matrix::print_values() {
	if (n > 10) {
		int fd_out = open("foutput.txt", O_RDWR | O_CREAT | O_TRUNC, 0777);
		dup2(fd_out, 1);
	}

	cout << "Values:" << endl;
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			cout.width(10);
			cout << values[i][j] << "   ";
		}
		cout << endl << endl;
	}
}

int main (int argc, char ** argv)
{
	if (argc < 4) {
		perror("No enough arguements");
		return -1;
	}

	/*Command line format: ./bin_file key [input_file_name] [n] mode_parameter. To fill the matrix via analytic function use the key "a" and also enter "n" parameter. To scan the matrix from the file use the key "f"  and also enter input_file_name.*/
	int n = 0;
	bool analytic = false;

	if (*argv[1] == 'a') {
		analytic = true;
		n = atoi(argv[2]);
	}

	else if (* argv[1] == 'f') {
		int fd = open (argv[2], O_RDWR, 0777);
        dup2 (fd, 0);
	}

	if (!analytic) cin >> n;

	Matrix A(n);

	if (analytic) A.set_values();
	else {
		A.get_values();
	}
	A.set_optim();
	int mode = atoi(argv[3]);

	A.eigenvalues(1e-8, mode);
	A.print_eigen_sollutions();
	A.print_values();
}
