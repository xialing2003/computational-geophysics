/*
Compile:
    h5c++ -I ../eigen heat_equation.cpp -o heat_equation
Run:
    ./heat_equation >> list.txt
This list.txt is a file to record the plotting steps. 
Use a brand-new file for each time interval choice.
*/

#include <bits/stdc++.h>
#include <../eigen/Eigen/Dense>
#include <../eigen/Eigen/SparseCholesky>
#include "H5Cpp.h"

using namespace Eigen;
using namespace H5;
using namespace std;

#define dx 1
#define D 1
double dt = 0.6 * dx * dx / D, t_record;
const int Nx = 101, Nt = int(100 / dt);

char output_1[] = "T_forward.h5";
char output_2[] = "T_backward.h5";
char output_3[] = "T_CNscheme.h5";

void output_slice(const VectorXd &T, const string &filename, double time) {
    H5File file(filename, H5F_ACC_RDWR);
    hsize_t dim[2] = {Nx, 1};
    DataSpace dataspace(2, dim);
    stringstream ss;
    ss << fixed << setprecision(2) << time;
    DataSet dataset = file.createDataSet("Time_" + ss.str(), PredType::NATIVE_DOUBLE, dataspace);
    dataset.write(T.data(), PredType::NATIVE_DOUBLE);
}

int main() {
    
    Matrix<double, Nx, 1> T;
    Eigen::SparseMatrix<double> A(Nx, Nx), B(Nx, Nx), C(Nx, Nx), E(Nx, Nx);

    for (int i = 1; i < Nx-1; i++) {
        A.insert(i, i - 1) = D * dt / (dx * dx);
        A.insert(i, i) = 1 - 2 * D * dt / (dx * dx);
        A.insert(i, i + 1) = D * dt / (dx * dx);
    }
    A.makeCompressed();

    B.insert(0,0) = 1;
    for (int i = 1; i < Nx-1; i++) {
        B.insert(i, i - 1) = - D * dt / (dx * dx);
        B.insert(i, i) = 1 + 2 * D * dt / (dx * dx);
        B.insert(i, i + 1) = - D * dt / (dx * dx);
    }
    B.insert(Nx-1, Nx-1) = 1;
    B.makeCompressed();

    C.insert(0,0) = 1;
    for (int i = 1; i < Nx-1; i++) {
        C.insert(i, i - 1) = - 0.5 * D * dt / (dx * dx);
        C.insert(i, i) = 1 + D * dt / (dx * dx);
        C.insert(i, i + 1) = - 0.5 * D * dt / (dx * dx);
    }
    C.insert(Nx-1, Nx-1) = 1;
    C.makeCompressed();

    E.setZero();
    E.insert(0,0) = 1;
    for (int i = 1; i < Nx-1; i++) {
        E.insert(i, i - 1) = 0.5 * D * dt / (dx * dx);
        E.insert(i, i) = 1 - D * dt / (dx * dx);
        E.insert(i, i + 1) = 0.5 * D * dt / (dx * dx);
    }
    E.insert(Nx-1, Nx-1) = 1;
    E.makeCompressed();

    // use the explicit scheme to solve this problem(forward difference in time)
    T.setZero();
    T(50) = 1;

    ofstream outFile_1(output_1, ios::out | ios::trunc);  // Open for writing, truncate existing content
    if (!outFile_1.is_open()) {
        cout << "Error opening the output_forward file!" << endl;
        return 0;
    }
    output_slice(T, output_1, 0);
    cout << "0" << endl;

    t_record = 10;
    for (int i = 1; i <= Nt; i++) {
        T = A * T;
        if(i *dt >= t_record){
            output_slice(T, output_1, i * dt);
            t_record += 10;
            cout << i << endl;  // using different schemes, I can still choose to plot at the same time step
        }
    }

    // use the implicit scheme to solve this problem(backward difference in time)
    T.setZero();
    T(50) = 1;

    ofstream outFile_2(output_2, ios::out | ios::trunc);  // Open for writing, truncate existing content
    if (!outFile_2.is_open()) {
        cout << "Error opening the output_backward file!" << endl;
        return 0;
    }
    output_slice(T, output_2, 0);

    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver_b;
    solver_b.compute(B);
    t_record = 10;
    for (int i = 1; i <= Nt; i++) {
        T = solver_b.solve(T);
        T(0) = 0;
        T(100) = 0;
        if(i *dt >= t_record){
            output_slice(T, output_2, i * dt);
            t_record += 10;
        }
    }

    // use the semi-implicit scheme to solve this problem (Crank-Nicolson scheme)
    T.setZero();
    T(50) = 1;

    ofstream outFile_3(output_3, ios::out | ios::trunc);  // Open for writing, truncate existing content
    if (!outFile_3.is_open()) {
        cout << "Error opening the output_backward file!" << endl;
        return 0;
    }
    output_slice(T, output_3, 0);

    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver_c;
    solver_c.compute(C);
    t_record = 10;
    for (int i = 1; i <= Nt; i++) {
        T = solver_c.solve(E * T);
        T(0) = 0;
        T(100) = 0;
        if(i *dt >= t_record){
            output_slice(T, output_3, i * dt);
            t_record += 10;
        }
    }

    return 0;
}
