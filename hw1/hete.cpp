#include <bits/stdc++.h>
#include <H5Cpp.h>
using namespace std;
const H5std_string FILE_NAME("hete_data.h5");

int main()
{
    int L = 100, rho = 1, c1 = 1, k1 = 1, c2 = 2, k2 = 4;
    double dx = 0.1, dt = 0.05;
    int si = L/dx, bi = 60/dx, sn = 2000;

    // Second Order System
    double u[sn + 1][si + 1], v2[sn + 1][si + 1]; //u_tn_xi
    for(int i = 0; i <= si; i++){
        u[0][i] = exp(-0.1*pow(dx*i-50, 2));
        u[1][i] = exp(-0.1*pow(dx*i-50, 2));
    }
    for(int n = 1; n < sn; n++){
        u[n+1][0] = 0;
        for(int i = 1; i <= bi; i++){
            u[n+1][i] = pow(c1*dt/dx, 2) * (u[n][i+1] - 2*u[n][i] + u[n][i-1]) + 2*u[n][i] - u[n-1][i];
        }
        for(int i = bi + 1; i <= si; i++){
            u[n+1][i] = pow(c2*dt/dx, 2) * (u[n][i+1] - 2*u[n][i] + u[n][i-1]) + 2*u[n][i] - u[n-1][i];
        }
        u[n+1][si] = u[n+1][si - 1];
    }

    for(int i = 0; i <= si; i++){
        v2[0][i] = 0;
        for(int n = 1; n < sn; n++){
            v2[n][i] = (u[n+1][i] - u[n-1][i]) / (2 * dt);
        }
        v2[sn][i] = (u[sn][i] - u[sn - 1][i]) / dt;
    }

    //improved First Order System
    double u_im[sn + 1][si + 1], v2_im[sn + 1][si + 1]; //u_tn_xi
    for(int i = 0; i <= si; i++){
        u_im[0][i] = exp(-0.1*pow(dx*i-50, 2));
        u_im[1][i] = exp(-0.1*pow(dx*i-50, 2));
    }
    for(int n = 1; n < sn; n++){
        u_im[n+1][0] = 0;
        for(int i = 1; i < bi; i++){
            u_im[n+1][i] = pow(c1*dt/dx, 2) * (u_im[n][i+1] - 2*u_im[n][i] + u_im[n][i-1]) + 2*u_im[n][i] - u_im[n-1][i];
        }
        u_im[n+1][bi] = pow(c2*dt/dx, 2) * (u_im[n][bi+1] - 2*u_im[n][bi] + u_im[n][bi-1]) + \
                        pow(dt/dx, 2) * (k2 - k1)/(2*rho) * (u_im[n][bi+1] - u_im[n][bi-1]) + \
                        2*u_im[n][bi] - u_im[n-1][bi];
        for(int i = bi + 1; i <= si; i++){
            u_im[n+1][i] = pow(c2*dt/dx, 2) * (u_im[n][i+1] - 2*u_im[n][i] + u_im[n][i-1]) + 2*u_im[n][i] - u_im[n-1][i];
        }
        u_im[n+1][si] = u_im[n+1][si - 1];
    }

    for(int i = 0; i <= si; i++){
        v2_im[0][i] = 0;
        for(int n = 1; n < sn; n++){
            v2_im[n][i] = (u_im[n+1][i] - u_im[n-1][i]) / (2 * dt);
        }
        v2_im[sn][i] = (u_im[sn][i] - u_im[sn - 1][i]) / dt;
    }

    // First Order System
    double u0[si + 1], v1[sn + 1][si + 1], T[sn + 1][si + 1];
    for(int i = 0; i <= si; i++){
        u0[i] = exp(-0.1*pow(dx*i-50, 2));
        v1[0][i] = 0;
    }
    for(int i = 1; i <= bi; i++){
        T[0][i] = k1 * (u0[i+1] - u0[i-1])/(2 * dx);
    }
    for(int i = bi + 1; i < si; i++){
        T[0][i] = k2 * (u0[i+1] - u0[i-1])/(2 * dx);
    }
    T[0][0] = k1 * (u0[1] - 0)/dx;
    T[0][si] = 0;

    for(int i = 1; i <= bi; i++){
        T[1][i] = (k1*dt/(2*dx)) * (v1[0][i+1] - v1[0][i-1]) + T[0][i];
        v1[1][i] = dt/(2*dx*rho) * (T[0][i+1] - T[0][i-1]) + v1[0][i];
    }
    for(int i = bi + 1; i < si; i++){
        T[1][i] = (k2*dt/(2*dx)) * (v1[0][i+1] - v1[0][i-1]) + T[0][i];
        v1[1][i] = dt/(2*dx*rho) * (T[0][i+1] - T[0][i-1]) + v1[0][i];
    }
    v1[1][0] = 0;
    T[1][0] = (k1*dt/dx) * (v1[1][1] - v1[1][0]) + T[0][0];
    T[1][si] = 0;
    v1[1][si] = dt/(dx*rho) * (T[0][si] - T[0][si-1]) + v1[0][si];

    for(int n = 1; n < sn; n++){
        v1[n+1][0] = 0;
        T[n+1][0] = (2*k1*dt/dx) * (v1[n][1] - v1[n][0]) + T[n-1][0];
        for(int i = 1; i <= bi; i++){
            T[n+1][i] = (k1*dt/dx) * (v1[n][i+1] - v1[n][i-1]) + T[n-1][i];
            v1[n+1][i] = dt/(dx*rho) * (T[n][i+1] - T[n][i-1]) + v1[n-1][i];
        }
        for(int i = bi + 1; i < si; i++){
            T[n+1][i] = (k2*dt/dx) * (v1[n][i+1] - v1[n][i-1]) + T[n-1][i];
            v1[n+1][i] = dt/(dx*rho) * (T[n][i+1] - T[n][i-1]) + v1[n-1][i];
        }
        T[n+1][si] = 0;
        v1[n+1][si] = 2*dt/(dx*rho) * (T[n][si] - T[n][si - 1]) + v1[n-1][si];
    }

    H5::H5File file(FILE_NAME, H5F_ACC_TRUNC);
    hsize_t dims1[2] = {hsize_t(sn + 1), hsize_t(si + 1)}; // data1 维度
    H5::DataSpace dataspace1(2, dims1);

    const char* dataset_names[] = {
        "2order_u", "2order_v", "1order_v","1order_T", 
        "2order_u_improve", "2order_v_improve"
    };

    double* dataset_data[] = {
        &u[0][0], &v2[0][0], &v1[0][0], &T[0][0], &u_im[0][0], &v2_im[0][0]
    };

    for (int i = 0; i < 6; ++i) {
        H5::DataSet dataset = file.createDataSet(dataset_names[i], H5::PredType::NATIVE_DOUBLE, dataspace1);
        dataset.write(dataset_data[i], H5::PredType::NATIVE_DOUBLE);
    }

    std::cout << "Multiple datasets written to " << FILE_NAME << " successfully!" << std::endl;

    return 0;
}