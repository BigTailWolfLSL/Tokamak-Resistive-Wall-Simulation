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
#include "RBmagnet.h"
#include "computeTimeStep.h"
#include "bili_interpolation.h"

using namespace std;

void getFlux(vector<vector<Equ>> &flux, vector<vector<Equ>> &u, double dx,
             double dt, char direction = 'x',int NUM_GROW=1) {
// 'dx' is already the delta in the 'direction'.
    int nx = flux.size();
    int ny = flux[0].size();
    if (direction == 'x') {
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny-1; j++) {
//                flux[i][j] = HLL_flux(u[NUM_GROW+i-1][NUM_GROW+j],u[NUM_GROW+i][NUM_GROW+j], dx, dt,direction);
                flux[i][j]= MUSCL_Hancock_flux(u,NUM_GROW+i,NUM_GROW+j,dx, dt,direction);
//TODO:We change the x-direction numerical method here.
            }
        }
    }
    if (direction == 'y') {
        for (int i = 0; i < nx-1; i++) {
            for (int j = 0; j < ny; j++) {
//                flux[i][j] = HLL_flux(u[NUM_GROW+i][NUM_GROW+j-1],u[NUM_GROW+i][NUM_GROW+j], dx, dt,direction);
                flux[i][j]= MUSCL_Hancock_flux(u,NUM_GROW+i,NUM_GROW+j,dx, dt,direction);
//TODO:We change the y-direction numerical method here.
            }
        }
    }
    return;
}

int main() {
    double x0 = 0.15;
    double x1 = 1.15;
    double y0 = 0.0;
    double y1 = 1.0;
    int x_nCells = 300;
    int y_nCells = 300;
    double tStart = 0.0;
    double tStop = 0.10;
    int NUM_GROW=2;//number of ghost cells
    SimMesh.set(x0,x1,x_nCells,y0,y1,y_nCells,NUM_GROW);
    double eta_wall=0.0342;

    double dx = SimMesh.dx;
    double dy = SimMesh.dy;
    double dt;
    vector<vector<Equ>> u(x_nCells + 2*NUM_GROW,
                          vector<Equ>(y_nCells + 2*NUM_GROW)); // 有效范围1~nCells+1
    vector<vector<Equ>> uRigidBody(x_nCells + 2*NUM_GROW,
                                   vector<Equ>(y_nCells + 2*NUM_GROW));
    vector<vector<Equ>> uPlus1(x_nCells + 2*NUM_GROW,
                               vector<Equ>(y_nCells + 2*NUM_GROW)); // 两个u之间有一个flux
    rigidBodyPhi(uPlus1,uRigidBody);
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

            u[i][j].rho = (x<=0.5?1.0:0.125);
            u[i][j].momentum_x=0;
            u[i][j].momentum_y=0;
            u[i][j].momentum_z = 0;
            u[i][j].pressure = (x<=0.5?1.0:0.1);
            u[i][j].B_x=0.75;
            u[i][j].B_y=0;
            u[i][j].B_z=(x<=0.5?1.0:-1.0);
            u[i][j].energy =u[i][j].pressure/(C_gamma-1.0)
                             +0.5*(u[i][j].momentum_x*u[i][j].momentum_x+u[i][j].momentum_y*u[i][j].momentum_y+u[i][j].momentum_z*u[i][j].momentum_z)/u[i][j].rho
                             +0.5*(u[i][j].B_x*u[i][j].B_x+u[i][j].B_y*u[i][j].B_y+u[i][j].B_z*u[i][j].B_z);

            uRigidBody[i][j].rho=0;
            uRigidBody[i][j].momentum_x = 0;
            uRigidBody[i][j].momentum_y = 0;
            uRigidBody[i][j].momentum_z = 0;
            uRigidBody[i][j].pressure = 0;
            uRigidBody[i][j].B_x=u[i][j].B_x;
            uRigidBody[i][j].B_y=u[i][j].B_y;
            uRigidBody[i][j].B_z=u[i][j].B_z;
            uRigidBody[i][j].energy =0;
        }
    }
    uPlus1=u;
    rigidBodyBC(u);
    transmissive(u,NUM_GROW);
    transmissive(uRigidBody,NUM_GROW);
    double t = tStart;
// 初始化结束

    while (t < tStop) {
        dt = ComputeTimeStep(u, dx, dy, 0.8,NUM_GROW);
        if (t + dt > tStop) {
            dt = tStop - t;
            t = tStop;
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
        u = uPlus1;
        rigidBodyBC(u);
        transmissive(u,NUM_GROW);
        //end of y update

        int NUM;NUM=ceil(eta_wall*20);if (NUM<1) NUM=1;
        for (int _=0;_<NUM;_++){
            setBC(uRigidBody, u);
            updateRigidBodyMagneticField(uRigidBody, dt / 2/NUM, eta_wall);
            transmissive(uRigidBody, NUM_GROW);
        }

        //x update
        getFlux(flux, u, dx, dt, 'x',NUM_GROW);
        for (int i = 0; i < x_nCells; i++) {
            for (int j = 0; j < y_nCells; j++) {
                uPlus1[NUM_GROW+i][NUM_GROW+j] = u[NUM_GROW+i][NUM_GROW+j] - (dt / dx) * (flux[i+1][j] - flux[i][j]);
            }
        }
        u = uPlus1;
        rigidBodyBC(u);
        transmissive(u,NUM_GROW);
        //end of x update

        for (int _=0;_<NUM;_++){
            setBC(uRigidBody, u);
            updateRigidBodyMagneticField(uRigidBody, dt / 2/NUM, eta_wall);
            transmissive(uRigidBody, NUM_GROW);
        }

        //beginning of y update
        getFlux(flux, u, dy, dt/2,'y',NUM_GROW);

        for (int i = 0; i < x_nCells; i++) {
            for (int j = 0; j < y_nCells; j++) {
                uPlus1[NUM_GROW+i][NUM_GROW+j] = u[NUM_GROW+i][NUM_GROW+j] - (dt/2 / dy) * (flux[i][j+1] - flux[i][j]);
            }
        }
        u = uPlus1;
        rigidBodyBC(u);
        transmissive(u,NUM_GROW);
        //end of y update

        cout<<"time ="<<t<<endl;

    }

    {
        ofstream outfile("/CFD_file/shockResistiveCylinder/Euler.dat");
        if (!outfile)
            return -1;

        for (int i = NUM_GROW; i < x_nCells + NUM_GROW; i++) {
            for (int j = NUM_GROW; j < y_nCells + NUM_GROW; j++) {
                double x = x0 + (i + 0.5 - NUM_GROW) * dx;
                double y = y0 + (j + 0.5 - NUM_GROW) * dy;
                if (get<0>(Phi[i][j])<=0) {
                    outfile << x << " " << y << " " << uRigidBody[i][j]
                            << " " << get<0>(Phi[i][j])
                            << endl;
                    continue;
                }
                outfile << x << " " << y << " " << u[i][j]
                        <<" "<<get<0>(Phi[i][j])
                        << endl;
            }
            outfile<<endl;
        }

        outfile.close();
    }

    return 0;
}
