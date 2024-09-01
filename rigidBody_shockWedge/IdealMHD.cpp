#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include <fstream>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <stdio.h>

#include "boundary.h"
#include "MHDvec.h"
#include "solver.h"
#include "mesh.h"
#include "rigidBodyBC.h"
#include "bili_interpolation.h"

using namespace std;

double ComputeTimeStep(vector<vector<Equ>> &u, double dx, double dy,
                       double C_CFL = 1,int NUM_GROW=1) {
    double a = 0;
    int nx = u.size();
    int ny = u[nx - 1].size();
    for (int i = NUM_GROW; i < nx - NUM_GROW; i++) {
        for (int j = NUM_GROW; j < ny - NUM_GROW; j++) {
            u[i][j].EoS();
            double at=max(abs(u[i][j].vel_x)+sqrt(0.5*(u[i][j].csSqu+u[i][j].caSqu+ sqrt((u[i][j].csSqu+u[i][j].caSqu)*(u[i][j].csSqu+u[i][j].caSqu)-4*u[i][j].csSqu*u[i][j].B_x*u[i][j].B_x/u[i][j].rho))),
                          abs(u[i][j].vel_y)+sqrt(0.5*(u[i][j].csSqu+u[i][j].caSqu+ sqrt((u[i][j].csSqu+u[i][j].caSqu)*(u[i][j].csSqu+u[i][j].caSqu)-4*u[i][j].csSqu*u[i][j].B_y*u[i][j].B_y/u[i][j].rho))));
            a=max(a,at);//TODO:change to 1D revise 2 place
        }
    }
    return C_CFL*min(dx,dy)/a;
}

void getFlux(vector<vector<Equ>> &flux, vector<vector<Equ>> &u, double dx,
             double dt, char direction = 'x',int NUM_GROW=1) {
// 'dx' is already the delta in the 'direction'.
    int nx = flux.size();
    int ny = flux[0].size();
    if (direction == 'x') {
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny-1; j++) {
//                flux[i][j] = HLLC_flux2(u[NUM_GROW+i-1][NUM_GROW+j],u[NUM_GROW+i][NUM_GROW+j], dx, dt,direction);
                flux[i][j]= MUSCL_Hancock_flux(u,NUM_GROW+i,NUM_GROW+j,dx, dt,direction);
//TODO:We change the x-direction numerical method here.
            }
        }
    }
    if (direction == 'y') {
        for (int i = 0; i < nx-1; i++) {
            for (int j = 0; j < ny; j++) {
//                flux[i][j] = HLLC_flux2(u[NUM_GROW+i][NUM_GROW+j-1],u[NUM_GROW+i][NUM_GROW+j], dx, dt,direction);
                flux[i][j]= MUSCL_Hancock_flux(u,NUM_GROW+i,NUM_GROW+j,dx, dt,direction);
//TODO:We change the y-direction numerical method here.
            }
        }
    }
    return;
}

double shockLocus(vector<vector<Equ>>& u, int i, int j){
    double T;
    double divergence=0;
    T=(u[i+1][j].rho-u[i-1][j].rho)/SimMesh.dx/2;
    divergence+=T*T;
    T=(u[i][j+1].rho-u[i][j-1].rho)/SimMesh.dy/2;
    divergence+=T*T;
    if (u[i][j].rho!=0) divergence=sqrt(divergence)/u[i][j].rho; else if (divergence!=0) divergence=1e10; else divergence=0;
    double adjustCoefficent=-20.0 /1000.0;
    double ans=exp(adjustCoefficent*divergence);
    return divergence;
}

int main() {
    double x0 = 0.0;
    double x1 = 23;
    double y0 = -5.7;
    double y1 = 5.7;
    int x_nCells = 920  *2;
    int y_nCells = 228*2  *2;
    double tStart = 0.0;
    double tStop = 6.9;

//    double x0 = 0.0;
//    double x1 = 15;
//    double y0 = -5.7;
//    double y1 = 5.7;
//    int x_nCells = 600 /2;
//    int y_nCells = 228*2 /2;
//    double tStart = 0.0;
//    double tStop = 1.9;

//double x0 = 0.0;
//double x1 = 1;
//double y0 = 0;
//double y1 = 1;
//int x_nCells = 128;
//int y_nCells = 128;
//double tStart = 0.0;
//double tStop = 0.25;

    int NUM_GROW=2;//number of ghost cells
    SimMesh.set(x0,x1,x_nCells,y0,y1,y_nCells,NUM_GROW);

    double dx = SimMesh.dx;
    double dy = SimMesh.dy;
    double dt;
    vector<vector<Equ>> u(x_nCells + 2*NUM_GROW,
                          vector<Equ>(y_nCells + 2*NUM_GROW)); // 有效范围1~nCells+1
    vector<vector<Equ>> uPlus1(x_nCells + 2*NUM_GROW,
                               vector<Equ>(y_nCells + 2*NUM_GROW)); // 两个u之间有一个flux
    rigidBodyPhi(uPlus1);
    vector<vector<Equ>> flux(x_nCells + 1,
                             vector<Equ>(y_nCells + 1)); // flux[i]=f_i+1/2
//    vector<vector<double>> divergence(x_nCells,vector<double>(y_nCells));
//    vector<vector<double>> psi(x_nCells+ 2*NUM_GROW,vector<double>(y_nCells+ 2*NUM_GROW));
//    vector<vector<vector<double>>> psi_flux(x_nCells+1,vector<vector<double>>(y_nCells+1,vector<double>(3)));
// 2D时，flux计算x方向上相当于x方向加了1/2，计算y方向上相当于y方向加了1/2
// 定义结束

    for (int i = NUM_GROW; i < x_nCells + NUM_GROW; i++) {
        for (int j = NUM_GROW; j < y_nCells + NUM_GROW; j++) {
            double x = x0 + (i + 0.5-NUM_GROW) * dx;
            double y = y0 + (j + 0.5-NUM_GROW) * dy;
//            u[i][j].rho = ((x+y)<800.0?1.0:0.125);
//            u[i][j].momentum_x = 0;
//            u[i][j].momentum_y = 0;
//            u[i][j].momentum_z = 0;
//            u[i][j].energy = ((x+y)<800.0?1.0:0.1)+0.78125;
//            u[i][j].B_x=0.53033-0.70710*((x+y)<800?1.0:-1.0);
//            u[i][j].B_y=0.53033+0.70710*((x+y)<800?1.0:-1.0);
//            u[i][j].B_z=0;

//            u[i][j].rho = (y<400?1.0:0.125);
//            u[i][j].momentum_x = 0;
//            u[i][j].momentum_y = 0;
//            u[i][j].momentum_z = 0;
//            u[i][j].energy = (y<400?1.0:0.1)+0.78125;
//            u[i][j].B_x=(y<400?-1.0:1.0);
//            u[i][j].B_y=0.75;
//            u[i][j].B_z=0;

//            double alpha=0;
//            double a_input=1.0,b_input=0,c_input=-0.5;
//            u[i][j].rho = (a_input*x+b_input*y+c_input<0?1.0:0.125);
//            u[i][j].momentum_x = 0;
//            u[i][j].momentum_y = 0;
//            u[i][j].momentum_z = 0;
//            u[i][j].pressure=(a_input*x+b_input*y+c_input<0?1.0:0.1);
////            double shockTube_Bx=0.75,shockTube_By=(a_input*x+b_input*y+c_input<0?1.0:-1.0);
//            double shockTube_Bx=0,shockTube_By=0;
//            u[i][j].B_x=shockTube_Bx*cos(alpha)+shockTube_By* cos(alpha+0.5*M_PI);
//            u[i][j].B_y=shockTube_Bx* sin(alpha)+shockTube_By* sin(alpha+0.5*M_PI);
//            u[i][j].B_z=0;
//            u[i][j].energy =u[i][j].pressure/(C_gamma-1)
//                             +0.5*(u[i][j].momentum_x*u[i][j].momentum_x+u[i][j].momentum_y*u[i][j].momentum_y+u[i][j].momentum_z*u[i][j].momentum_z)/u[i][j].rho
//                             +0.5*(u[i][j].B_x*u[i][j].B_x+u[i][j].B_y*u[i][j].B_y+u[i][j].B_z*u[i][j].B_z);

//            double M=1.3;
////            double rho_0=1.225e-3,p_0=1.013e-4;
//            double rho_0=1.225,p_0=101325;
//            double rho_s=rho_0* ((C_gamma+1)*M*M)/(2+(C_gamma-1)*M*M);
//            double p_s=p_0* (2*C_gamma*M*M-(C_gamma-1))/(C_gamma+1);
//            double u_s= sqrt(  (rho_s-rho_0)*(p_s-p_0)/(rho_s*rho_0)  );
//            u[i][j].rho = (x<9?rho_s:rho_0);
//            u[i][j].momentum_x = (x<9?u_s:0);
//            u[i][j].momentum_y = 0;
//            u[i][j].momentum_z = 0;
//            u[i][j].pressure=(x<9?p_s:p_0);
//            u[i][j].B_x=0;
//            u[i][j].B_y=0;
//            u[i][j].B_z=0;
//            u[i][j].energy =u[i][j].pressure/(C_gamma-1)
//                             +0.5*(u[i][j].momentum_x*u[i][j].momentum_x+u[i][j].momentum_y*u[i][j].momentum_y+u[i][j].momentum_z*u[i][j].momentum_z)/u[i][j].rho
//                             +0.5*(u[i][j].B_x*u[i][j].B_x+u[i][j].B_y*u[i][j].B_y+u[i][j].B_z*u[i][j].B_z);

            double M=1.3;
            //            double rho_0=1.225e-3,p_0=1.013e-4;
            double rho_0=1.225,p_0=1.01325;
            double rho_s=rho_0* ((C_gamma+1)*M*M)/(2+(C_gamma-1)*M*M);
            double p_s=p_0* (2*C_gamma*M*M-(C_gamma-1))/(C_gamma+1);
            double u_s= sqrt(  (rho_s-rho_0)*(p_s-p_0)/(rho_s*rho_0)  );
            u[i][j].rho = (x<9.3?rho_s:rho_0);
            u[i][j].momentum_x = (x<9.3?u_s:0);
            u[i][j].momentum_y = 0;
            u[i][j].momentum_z = 0;
            u[i][j].pressure=(x<9.3?p_s:p_0);
            u[i][j].B_x=0;
            u[i][j].B_y=0;
            u[i][j].B_z=0;
            u[i][j].energy =u[i][j].pressure/(C_gamma-1)
                             +0.5*(u[i][j].momentum_x*u[i][j].momentum_x+u[i][j].momentum_y*u[i][j].momentum_y+u[i][j].momentum_z*u[i][j].momentum_z)/u[i][j].rho
                             +0.5*(u[i][j].B_x*u[i][j].B_x+u[i][j].B_y*u[i][j].B_y+u[i][j].B_z*u[i][j].B_z);
        }
    }
    uPlus1=u;
    rigidBodyBC(u);
//    diagonal_transmissive(u,NUM_GROW);
    customBoundary3(u,NUM_GROW);

    double t = tStart;

    double checkpoint[]={2.70,3.4,4.4,4.9,5.4,5.7,5.80,6.6,6.7,6.8,9.1,100001};
    int l=0;
// 初始化结束

    while (t < tStop) {
        dt = ComputeTimeStep(u, dx, dy, 0.4,NUM_GROW);
        if (t + dt > tStop) {
            dt = tStop - t;
            t = tStop;
        } else if (t<0.5&&(t+dt)>=0.5){
            t = t + dt;
            double thre=9.72222222;
            Equ T=u[(int)(thre/23.0*x_nCells) + NUM_GROW][y_nCells/2+NUM_GROW];
            for (int i = NUM_GROW; i < (int)(thre/23.0*x_nCells) + NUM_GROW; i++) {
                for (int j = NUM_GROW; j < y_nCells + NUM_GROW; j++) {
                    u[i][j]=T;
                }
            }
            customBoundary3(u,NUM_GROW);
        } else if (t<checkpoint[l]&&(t+dt)>=checkpoint[l]){
            dt=checkpoint[l]-t;
            t=checkpoint[l];
        } else {
            t = t + dt;
        }


        //beginning of y update
        getFlux(flux, u, dy, dt/2,'y',NUM_GROW);
        for (int i = 0; i < x_nCells; i++) {
            for (int j = 0; j < y_nCells; j++) {
                uPlus1[NUM_GROW+i][NUM_GROW+j] = u[NUM_GROW+i][NUM_GROW+j] - (dt/2 / dy) * (flux[i][j+1] - flux[i][j]);
            }
        }
//        diagonal_transmissive(uPlus1,NUM_GROW);
        u = uPlus1;
        rigidBodyBC(u);
        customBoundary3(u,NUM_GROW);
        //end of y update


        //x update
        getFlux(flux, u, dx, dt, 'x',NUM_GROW);
        for (int i = 0; i < x_nCells; i++) {
            for (int j = 0; j < y_nCells; j++) {
                uPlus1[NUM_GROW+i][NUM_GROW+j] = u[NUM_GROW+i][NUM_GROW+j] - (dt / dx) * (flux[i+1][j] - flux[i][j]);
            }
        }
//        diagonal_transmissive(uPlus1,NUM_GROW);
        u = uPlus1;
        rigidBodyBC(u);
        customBoundary3(u,NUM_GROW);
        //end of x update


        //beginning of y update
        getFlux(flux, u, dy, dt/2,'y',NUM_GROW);

        for (int i = 0; i < x_nCells; i++) {
            for (int j = 0; j < y_nCells; j++) {
                uPlus1[NUM_GROW+i][NUM_GROW+j] = u[NUM_GROW+i][NUM_GROW+j] - (dt/2 / dy) * (flux[i][j+1] - flux[i][j]);
            }
        }
//        diagonal_transmissive(uPlus1,NUM_GROW);
        u = uPlus1;
        rigidBodyBC(u);
        customBoundary3(u,NUM_GROW);
        //end of y update

        cout<<"time ="<<t<<endl;

        if (t==checkpoint[l]){
            string out1="\\\\wsl.localhost\\Ubuntu-18.04\\CFD_file\\rigidBody\\shock_wedge_";
            string out2=to_string((int)(checkpoint[l]*100));
            string out3=".dat";
            string outf=out1+out2+out3;
            ofstream outfile(outf);
            if (!outfile)
                return -1;

            for (int i = NUM_GROW; i < x_nCells + NUM_GROW; i++) {
                for (int j = NUM_GROW; j < y_nCells + NUM_GROW; j++) {
                    double x = x0 + (i + 0.5 - NUM_GROW) * dx;
                    double y = y0 + (j + 0.5 - NUM_GROW) * dy;
                    if (get<0>(Phi[i][j])<0) {
                        Equ t(1.0e+5, 1.0e+10, 1.0e+10, 1.0e+10, 1.0e+15, 1.0e+5, 1.0e+5, 1.0e+5);
                        outfile << x << " " << y << " " << t
                                <<get<0>(Phi[i][j])<<" "
                                << 100
                                << endl;
                        continue;
                    }
                    outfile << x << " " << y << " " << u[i][j]
                            <<get<0>(Phi[i][j])<<" "
                            << shockLocus(u,i,j)
                            << endl;
                }
                outfile<<endl;
            }

            outfile.close();
            if (checkpoint[l]>100000) break; else l++;
        }
    }

    {
        ofstream outfile("\\\\wsl.localhost\\Ubuntu-18.04\\CFD_file\\rigidBody\\shock_wedge.dat");
        if (!outfile)
            return -1;

        for (int i = NUM_GROW; i < x_nCells + NUM_GROW; i++) {
            for (int j = NUM_GROW; j < y_nCells + NUM_GROW; j++) {
                double x = x0 + (i + 0.5 - NUM_GROW) * dx;
                double y = y0 + (j + 0.5 - NUM_GROW) * dy;
                if (get<0>(Phi[i][j])<0) {
                    Equ t(1.0e+5, 1.0e+10, 1.0e+10, 1.0e+10, 1.0e+15, 1.0e+5, 1.0e+5, 1.0e+5);
                    outfile << x << " " << y << " " << t
                            <<get<0>(Phi[i][j])<<" "
                            << 40
                            << endl;
                    continue;
                }
                outfile << x << " " << y << " " << u[i][j]
                        <<get<0>(Phi[i][j])<<" "
                        << shockLocus(u,i,j)
                        << endl;
            }
            outfile<<endl;
        }

        outfile.close();
    }

//    {
//        ofstream outfile2("\\\\wsl.localhost\\Ubuntu-18.04\\CFD_file\\rigidBody\\Euler_lineout.dat");
//        if (!outfile2)
//            return -1;
//
//        int resolution = 1000;
//        double output_dx = (outputEnd_x - outputBegin_x) / (double) resolution;
//        double output_dy = (outputEnd_y - outputBegin_y) / (double) resolution;
//        for (int i = 0; i < resolution; i++) {
//            double x = outputBegin_x + i * output_dx;
//            double y = outputBegin_y + i * output_dy;
//            double cell_x = (x - x0) / dx - 0.5 + NUM_GROW;
//            double cell_y = (y - y0) / dy - 0.5 + NUM_GROW;
//            int x_1 = static_cast<int>(cell_x), y_1 = static_cast<int>(cell_y);
//            int x_2 = x_1 + 1, y_2 = y_1 + 1;
//            Equ result = bilinearInterpolation(u[x_1][y_1], u[x_2][y_1], u[x_1][y_2], u[x_2][y_2], (cell_x - x_1), (cell_y - y_1));
//            outfile2 << (1.5 * i / resolution) << " " << result
//                     << endl;
//        }
//        outfile2.close();
//    }

    return 0;
}
