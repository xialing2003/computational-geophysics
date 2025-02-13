#include <bits/stdc++.h>
using namespace std;
int main()
{
    double u[501][1005]; //u_tn_xi
    int c = 1, rho = 1, k = 1, L = 100;
    double dx = 0.1, dt = 0.1;

    //initialization
    for(int i = 0; i<=1000; i++){
        u[0][i] = exp(-0.1*pow(dx*i-50, 2));
    }
    for(int n = 2; n <= 500; n++){
        for(int i = 0; i <= 1000; i++){
            u[n+1][i] = pow(c*dt/dx, 2) * (u[n][i+1] - 2*u[n][i] + u[n][i-1]) + 2*u[n][i] - u[n-1][i];
        }
    }
    return 0;
}