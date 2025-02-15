#include <bits/stdc++.h>
#include <H5Cpp.h>
using namespace std;
const H5std_string FILE_NAME("homo_data.h5");

int main()
{
    int L = 100, rho = 1, c = 1, k = 1;
    int si = 1000, sn = 1000;
    double dx = 0.1, dt = 0.1;

    // Second Order System
    double u_d[sn + 1][si + 1], u_n[sn + 1][si + 1]; //u_tn_xi
    for(int i = 0; i <= si; i++){
        u_d[0][i] = exp(-0.1*pow(dx*i-50, 2));
        u_d[1][i] = exp(-0.1*pow(dx*i-50, 2));
        u_n[0][i] = exp(-0.1*pow(dx*i-50, 2));
        u_n[1][i] = exp(-0.1*pow(dx*i-50, 2));
    }
    //// Dirichlet boundary conditions u(0,t) = 0
    for(int n = 1; n < sn; n++){
        u_d[n+1][0] = 0;
        u_d[n+1][si] = 0;
        for(int i = 1; i < si; i++){
            u_d[n+1][i] = pow(c*dt/dx, 2) * (u_d[n][i+1] - 2*u_d[n][i] + u_d[n][i-1]) + 2*u_d[n][i] - u_d[n-1][i];
        }
    }
    //// Neumann boundary conditions ux(0,t) = 0
    for(int n = 1; n < sn; n++){
        for(int i = 1; i < si; i++){
            u_n[n+1][i] = pow(c*dt/dx, 2) * (u_n[n][i+1] - 2*u_n[n][i] + u_n[n][i-1]) + 2*u_n[n][i] - u_n[n-1][i];
        }
        u_n[n+1][0] = u_n[n+1][1];
        u_n[n+1][si] = u_n[n+1][si - 1];
    }
    //// get corresponding velocity
    double v2_d[sn + 1][si + 1], v2_n[sn + 1][si + 1];

    for(int i = 0; i <= si; i++){
        v2_d[0][i] = 0;
        v2_n[0][i] = 0;
        for(int n = 1; n < sn; n++){
            v2_d[n][i] = (u_d[n+1][i] - u_d[n-1][i]) / (2 * dt);
            v2_n[n][i] = (u_n[n+1][i] - u_n[n-1][i]) / (2 * dt);
        }
        v2_d[sn][i] = (u_d[sn][i] - u_d[sn - 1][i]) / dt;
        v2_n[sn][i] = (u_n[sn][i] - u_n[sn - 1][i]) / dt;
    }

    // First Order System
    double u0[si + 1], v1_d[sn + 1][si + 1], T_d[sn + 1][si + 1], v1_n[sn + 1][si + 1], T_n[sn + 1][si + 1];
    for(int i = 0; i <= si; i++){
        u0[i] = exp(-0.1*pow(dx*i-50, 2));
        v1_d[0][i] = 0;
        v1_n[0][i] = 0;
    }
    //// Dirichlet boundary conditions u(0,t) = 0
    for(int i = 1; i < si; i++){
        T_d[0][i] = k * (u0[i+1] - u0[i-1])/(2 * dx);
    }
    T_d[0][0] = k * (u0[1] - 0)/dx;
    T_d[0][si] = k * (0 - u0[si - 1])/dx;

    for(int i = 1; i < si; i++){
        T_d[1][i] = (k*dt/(2*dx)) * (v1_d[0][i+1] - v1_d[0][i-1]) + T_d[0][i];
        v1_d[1][i] = dt/(2*dx*rho) * (T_d[0][i+1] - T_d[0][i-1]) + v1_d[0][i];
    }
    v1_d[1][0] = 0;
    v1_d[1][si] = 0;
    T_d[1][0] = (k*dt/dx) * (v1_d[1][1] - v1_d[1][0]) + T_d[0][0];
    T_d[1][si] = (k*dt/dx) * (v1_d[1][si] - v1_d[1][si - 1]) + T_d[1][si - 1];

    for(int n = 1; n < sn; n++){
        v1_d[n+1][0] = 0;
        v1_d[n+1][si] = 0;
        T_d[n+1][0] = (2*k*dt/dx) * (v1_d[n][1] - v1_d[n][0]) + T_d[n-1][0];
        T_d[n+1][si] = (2*k*dt/dx) * (v1_d[n][si] - v1_d[n][si - 1]) + T_d[n-1][si];
        for(int i = 1; i < si; i++){
            T_d[n+1][i] = (k*dt/dx) * (v1_d[n][i+1] - v1_d[n][i-1]) + T_d[n-1][i];
            v1_d[n+1][i] = dt/(dx*rho) * (T_d[n][i+1] - T_d[n][i-1]) + v1_d[n-1][i];
        }
    }

    //// Neumann boundary conditions ux(0,t) = 0
    for(int i = 1; i < si; i++){
        T_n[0][i] = k * (u0[i+1] - u0[i-1])/(2 * dx);
    }
    T_n[0][0] = 0;
    T_n[0][si] = 0;

    for(int i = 1; i < si; i++){
        T_n[1][i] = (k*dt/(2*dx)) * (v1_n[0][i+1] - v1_n[0][i-1]) + T_n[0][i];
        v1_n[1][i] = dt/(2*dx*rho) * (T_n[0][i+1] - T_n[0][i-1]) + v1_n[0][i];
    }
    T_n[1][0] = 0;
    T_n[1][si] = 0;
    v1_n[1][0] = dt/(dx*rho) * (T_n[0][1] - T_n[0][0]) + v1_n[0][0];
    v1_n[1][si] = dt/(dx*rho) * (T_n[0][si] - T_n[0][si-1]) + v1_n[0][si];
    
    for(int n = 1; n < sn; n++){
        T_n[n+1][0] = 0;
        T_n[n+1][si] = 0;
        v1_n[n+1][0] = 2*dt/(dx*rho) * (T_n[n][1] - T_n[n][0]) + v1_n[n-1][0];
        v1_n[n+1][si] = 2*dt/(dx*rho) * (T_n[n][si] - T_n[n][si - 1]) + v1_n[n-1][si];
        for(int i = 1; i < 1000; i++){
            T_n[n+1][i] = (k*dt/dx) * (v1_n[n][i+1] - v1_n[n][i-1]) + T_n[n-1][i];
            v1_n[n+1][i] = dt/(dx*rho) * (T_n[n][i+1] - T_n[n][i-1]) + v1_n[n-1][i];
        }
    }

    H5::H5File file(FILE_NAME, H5F_ACC_TRUNC);
    hsize_t dims1[2] = {hsize_t(sn + 1), hsize_t(si + 1)}; // data1 维度
    H5::DataSpace dataspace1(2, dims1);

    const char* dataset_names[] = {
        "2order_dirich_u", "2order_neumann_u",
        "2order_dirich_v", "2order_neumann_v",
        "1order_dirich_v", "1order_neumann_v",
        "1order_dirich_T", "1order_neumann_T"
    };

    double* dataset_data[] = {
        &u_d[0][0], &u_n[0][0], &v2_d[0][0], &v2_n[0][0],
        &v1_d[0][0], &v1_n[0][0], &T_d[0][0], &T_n[0][0]
    };

    for (int i = 0; i < 8; ++i) {
        H5::DataSet dataset = file.createDataSet(dataset_names[i], H5::PredType::NATIVE_DOUBLE, dataspace1);
        dataset.write(dataset_data[i], H5::PredType::NATIVE_DOUBLE);
    }

    std::cout << "Multiple datasets written to " << FILE_NAME << " successfully!" << std::endl;

    return 0;
}