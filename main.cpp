#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <chrono>
#include <sstream>
#include <string>

double** matrix;;

void
create_matrix(int& size_of_matrix)
{
    std::ifstream file;
    file.open("/home/l.burtsev/Desktop/CHMI/in.csv");
    std::string line;
    while(getline(file, line)) {
        size_of_matrix++;
    }

    file.close();
    file.open("/home/l.burtsev/Desktop/CHMI/in.csv");
    matrix = new double*[size_of_matrix];
    std::string x, num;
    for(int i = 0; i < size_of_matrix; i++) {
        matrix[i] = new double[size_of_matrix];
        getline(file, line);
        std::stringstream str(line);
        for(int j = 0; j < size_of_matrix; j++) {
            getline(str, num, ',');
            matrix[i][j] = std::stod(num);
        }
    }
    file.close();
}

void
create_true_solve(double* true_solve, int size_of_matrix)
{
    for(int i = 0; i < size_of_matrix; i++) {
        true_solve[i] = 2 * (double)rand() / RAND_MAX - 1;
    }
}

void
create_matrix_F(double** matrix, double* F, double* true_solve, int size_of_matrix)
{
    for(int i = 0; i < size_of_matrix; i++) {
        F[i] = 0;
        for(int j = 0; j < size_of_matrix; j++) {
            F[i] += matrix[i][j] * true_solve[j];
        }
    }
}

void
create_matrix_P(double** matrix, double** P, int size_of_matrix, int start_position)
{
    double* u = new double[size_of_matrix - start_position];
    double sum = 0;
    for(int i = start_position; i < size_of_matrix; i++) {
        u[i - start_position] = matrix[i][start_position];
        sum += u[i - start_position] * u[i - start_position];
    }
    u[0] -= sqrt(sum);
    sum -= matrix[start_position][start_position] * matrix[start_position][start_position] - u[0] * u[0];
    sum /= 2;

    for(int i = 0; i < size_of_matrix; i++) {
        for(int j = 0; j < size_of_matrix; j++) {
            if(i < start_position || j < start_position) {
                if(i == j) {
                    P[i][j] = 1;
                } else {
                    P[i][j] = 0;
                }
            } else {
                if(i == j) {
                    P[i][j] = 1 - u[i -start_position] * u[i - start_position] / sum;
                } else {
                    P[i][j] = -u[i - start_position] * u[j - start_position] / sum;
                }
            }
        }
    }
}

void
multiplicate_P_matrix(double** P, double** matrix, int size_of_matrix)
{
    double** new_matrix = new double*[size_of_matrix];

    for (int i = 0; i < size_of_matrix; i++) {
        new_matrix[i] = new double[size_of_matrix];
        for (int j = 0; j < size_of_matrix; j++) {
            new_matrix[i][j] = 0;
            for (int k = 0; k < size_of_matrix; k++) {
                new_matrix[i][j] += P[i][k] * matrix[k][j];
            }
        }
    }

    for (int i = 0; i < size_of_matrix; i++) {
        for (int j = 0; j < size_of_matrix; j++) {
            matrix[i][j] = new_matrix[i][j];
        }
    }
}

void
multiplicate_P_F(double** P, double* F, int size_of_matrix)
{
    double* new_vec = new double[size_of_matrix];

    for(int i = 0; i < size_of_matrix; i++) {
        new_vec[i] = 0;
        for(int j = 0; j < size_of_matrix; j++) {
            new_vec[i] += P[i][j] * F[j];
        }
    }

    for(int i = 0; i < size_of_matrix; i++) {
        F[i] = new_vec[i];
    }
}

int
main(void)
{
    int size_of_matrix = 0;
    create_matrix(size_of_matrix);
    double* true_solve = new double[size_of_matrix];
    double* F = new double[size_of_matrix];
    double** P = new double*[size_of_matrix];
    for(int i = 0; i < size_of_matrix; i++) {
        P[i] = new double[size_of_matrix];
    }
    double* solve = new double[size_of_matrix];
    create_true_solve(true_solve, size_of_matrix);
    create_matrix_F(matrix, F, true_solve, size_of_matrix);

    auto start = std::chrono::high_resolution_clock::now();
    for(int i = 0; i < size_of_matrix - 1; i++) {
        create_matrix_P(matrix, P, size_of_matrix, i);
        multiplicate_P_matrix(P, matrix, size_of_matrix);
        multiplicate_P_F(P, F, size_of_matrix);
    }

    double sum;
    for(int i = size_of_matrix - 1; i >= 0; i--) {
        sum = 0;
        for(int j = size_of_matrix - 1; j > i; j--) {
            sum += matrix[i][j] * solve[j];
        }
        solve[i] = (F[i] - sum) / matrix[i][i];
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Время работы программы: " << duration.count() << " миллисекунд" << std::endl;
    
    double pogr = 0;
    for(int i = 0; i < size_of_matrix; i++) {
        if(fabs(true_solve[i] - solve[i]) > pogr) {
            pogr = fabs(true_solve[i] - solve[i]);
        }
    }
    std::cout << "Максимум норма погрешности: " << pogr << std::endl;
    return 0;
}
