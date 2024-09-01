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
#include "bili_interpolation.h"

using namespace std;

double ComputeTimeStep(vector<vector<Equ>> &u, double dx, double dy,
                       double C_CFL = 1,int NUM_GROW=1) {
    double a = 0;
    int nx = u.size();
    int ny = u[nx - 1].size();
    for (int i = NUM_GROW; i < nx - NUM_GROW; i++) {
        for (int j = NUM_GROW; j < ny - NUM_GROW; j++) {
            if (get<0>(Phi[i][j])<0){ continue;}
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
    double x0 = 0.0;
    double x1 = 2.0;
    double y0 = 0.0;
    double y1 = 2.0;
    int x_nCells = 160;
    int y_nCells = 160;
    double tStart = 0.0;
    double tStop = 0.1;
    int NUM_GROW=2;//number of ghost cells
    SimMesh.set(x0,x1,x_nCells,y0,y1,y_nCells,NUM_GROW);
    double eta_wall=0.1;
    setEta(eta_wall);

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

//            double alpha=0.0;
//            double a_input=1.0,b_input=0.0,c_input=-1.0;
            u[i][j].rho = (a_input*x+b_input*y+c_input<0?1.0:0.125);
            u[i][j].momentum_x = 0;
            u[i][j].momentum_y = 0;
            u[i][j].momentum_z = 0;
            u[i][j].pressure = (a_input*x+b_input*y+c_input<0?1.0:0.1);
            double shockTube_Bx=0.75,shockTube_By=(a_input*x+b_input*y+c_input<0?1.0:-1.0);
//            double shockTube_Bx=0,shockTube_By=0;
            u[i][j].B_x=shockTube_Bx*cos(alpha);
            u[i][j].B_y=shockTube_Bx*sin(alpha);
            u[i][j].B_z=shockTube_By;
            u[i][j].energy =u[i][j].pressure/(C_gamma-1)
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

//        int NUM;NUM=ceil(eta_wall*20);if (NUM<1) NUM=1;
//        for (int _=0;_<NUM;_++){
//            setBC(uRigidBody, u);
//            updateRigidBodyMagneticField(uRigidBody, dt / 2/NUM, eta_wall);
//            transmissive(uRigidBody, NUM_GROW);
//        }
        static const double resistiveSubDt=0.8/eta_wall/2.0/(1/SimMesh.dx/SimMesh.dx+1/SimMesh.dy/SimMesh.dy);
        {
            double subt=0;
            double subDt;
            double subtStop=dt/2;
            while (subt<subtStop){
                if (subt+resistiveSubDt>subtStop){
                    subDt=subtStop-subt;
                    subt=subtStop;
                } else {
                    subDt=resistiveSubDt;
                    subt=subt+resistiveSubDt;
                }
                setBC(uRigidBody, u);
                updateRigidBodyMagneticField(uRigidBody, subDt, eta_wall);
                transmissive(uRigidBody, NUM_GROW);
            }
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

//        for (int _=0;_<NUM;_++){
//            setBC(uRigidBody, u);
//            updateRigidBodyMagneticField(uRigidBody, dt / 2/NUM, eta_wall);
//            transmissive(uRigidBody, NUM_GROW);
//        }
        {
            double subt=0;
            double subDt;
            double subtStop=dt/2;
            while (subt<subtStop){
                if (subt+resistiveSubDt>subtStop){
                    subDt=subtStop-subt;
                    subt=subtStop;
                } else {
                    subDt=resistiveSubDt;
                    subt=subt+resistiveSubDt;
                }
                setBC(uRigidBody, u);
                updateRigidBodyMagneticField(uRigidBody, subDt, eta_wall);
                transmissive(uRigidBody, NUM_GROW);
            }
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
        ofstream outfile("/CFD_file/ResistiveRotated_BrioWu/RotatedBrioWu.dat");
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

    {
        ofstream outfile2("/CFD_file/ResistiveRotated_BrioWu/RotatedBrioWu_lineout.dat");
        if (!outfile2)
            return -1;

        int resolution = 1000;
        double shifting_x=0;//(0.5-4*(dx*sin(alpha)+dy*cos(alpha)))* sin(alpha);
        double shifting_y=0;//(0.5-4*(dx*sin(alpha)+dy*cos(alpha)))* (-cos(alpha));
        double output_dx = (outputEnd_x - outputBegin_x) / (double) resolution;
        double output_dy = (outputEnd_y - outputBegin_y) / (double) resolution;
        for (int i = 0; i < resolution; i++) {
            double x = outputBegin_x + i * output_dx+shifting_x;
            double y = outputBegin_y + i * output_dy+shifting_y;
            double cell_x = (x - x0) / dx - 0.5 + NUM_GROW;
            double cell_y = (y - y0) / dy - 0.5 + NUM_GROW;
            int x_1 = static_cast<int>(cell_x), y_1 = static_cast<int>(cell_y);
            int x_2 = x_1 + 1, y_2 = y_1 + 1;
            Equ result = bilinearInterpolation(u[x_1][y_1], u[x_2][y_1], u[x_1][y_2], u[x_2][y_2], (cell_x - x_1), (cell_y - y_1));
            outfile2 << (1.5 * i / resolution) << " " << result
                     << endl;
        }
        outfile2.close();
    }
    {
        ofstream outfile2("/CFD_file/ResistiveRotated_BrioWu/RotatedBrioWu_lineout_boundary.dat");
        if (!outfile2)
            return -1;

        int resolution = 1000;
        double shifting_x=(0.5-1*(dx*sin(alpha)+dy*cos(alpha)))* sin(alpha);
        double shifting_y=(0.5-1*(dx*sin(alpha)+dy*cos(alpha)))* (-cos(alpha));
        double output_dx = (outputEnd_x - outputBegin_x) / (double) resolution;
        double output_dy = (outputEnd_y - outputBegin_y) / (double) resolution;
        for (int i = 0; i < resolution; i++) {
            double x = outputBegin_x + i * output_dx+shifting_x;
            double y = outputBegin_y + i * output_dy+shifting_y;
            double cell_x = (x - x0) / dx - 0.5 + NUM_GROW;
            double cell_y = (y - y0) / dy - 0.5 + NUM_GROW;
            int x_1 = static_cast<int>(cell_x), y_1 = static_cast<int>(cell_y);
            int x_2 = x_1 + 1, y_2 = y_1 + 1;
            Equ result = bilinearInterpolation(u[x_1][y_1], u[x_2][y_1], u[x_1][y_2], u[x_2][y_2], (cell_x - x_1), (cell_y - y_1));
            outfile2 << (1.5 * i / resolution) << " " << result
                     << endl;
        }
        outfile2.close();
    }

    return 0;
}
