#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include <fstream>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <stdio.h>
#include <omp.h>

#include "boundary.h"
#include "MHDvec.h"
#include "mesh.h"
#include "rigidBodyBC.h"
#include "bili_interpolation.h"
#include "solver.h"


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
            a=max(a,at);//TODO:DIM 2places to revise
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

int main() {
    double x0 = 0.0;
    double x1 = 2.0;
    double y0 = 0.0;
    double y1 = 2.0;
    int x_nCells = 100;
    int y_nCells = 100;
    double tStart = 0.0;
    double tStop = 0.20;
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


//    vector<vector<double>> divergence(x_nCells,vector<double>(y_nCells));//TODO:div
//    vector<vector<double>> psi(x_nCells+ 2*NUM_GROW,vector<double>(y_nCells+ 2*NUM_GROW));
//    vector<vector<vector<double>>> psi_flux(x_nCells+1,vector<vector<double>>(y_nCells+1,vector<double>(3)));


// 2D时，flux计算x方向上相当于x方向加了1/2，计算y方向上相当于y方向加了1/2
// 定义结束

    for (int i = NUM_GROW; i < x_nCells + NUM_GROW; i++) {
        for (int j = NUM_GROW; j < y_nCells + NUM_GROW; j++) {
            double x = SimMesh.getx(i);
            double y = SimMesh.gety(j);

//            u[i][j].rho = 1;
//            u[i][j].momentum_x = 1;
//            u[i][j].momentum_y = 0;
//            u[i][j].momentum_z = 0;
//            u[i][j].energy =1;
//            u[i][j].B_x=0;
//            u[i][j].B_y=0;
//            u[i][j].B_z=0;

//            u[i][j].rho = (x+y<2?1.0:0.125);
//            u[i][j].momentum_x = 0;
//            u[i][j].momentum_y = 0;
//            u[i][j].momentum_z = 0;
//            u[i][j].energy = (y+x<2?2.5:0.25);
//            u[i][j].B_x=0;
//            u[i][j].B_y=0;
//            u[i][j].B_z=0;

//            u[i][j].rho = (a_input*x+b_input*y+c_input>0?1.0:0.125);
//            u[i][j].momentum_x = 0;
//            u[i][j].momentum_y = 0;
//            u[i][j].momentum_z = 0;
//            u[i][j].energy = (a_input*x+b_input*y+c_input>0?2.5:0.25);
//            u[i][j].B_x=0;
//            u[i][j].B_y=0;
//            u[i][j].B_z=0;

//            u[i][j].rho = (a_input*x+b_input*y+c_input>0?1.0:0.125);
//            u[i][j].momentum_x = 0;
//            u[i][j].momentum_y = 0;
//            u[i][j].momentum_z = 0;
//            u[i][j].pressure=(a_input*x+b_input*y+c_input>0?1.0:0.1);
//            double shockTube_Bx=0.75,shockTube_By=(a_input*x+b_input*y+c_input>0?1.0:-1.0);
//            u[i][j].B_x=shockTube_Bx*cos(alpha)+shockTube_By* cos(alpha+0.5*M_PI);
//            u[i][j].B_y=shockTube_Bx* sin(alpha)+shockTube_By* sin(alpha+0.5*M_PI);
//            u[i][j].B_z=0;
//            u[i][j].energy = u[i][j].pressure/(C_gamma-1)
//                             +0.5*(u[i][j].momentum_x*u[i][j].momentum_x+u[i][j].momentum_y*u[i][j].momentum_y+u[i][j].momentum_z*u[i][j].momentum_z)/u[i][j].rho
//                             +0.5*(u[i][j].B_x*u[i][j].B_x+u[i][j].B_y*u[i][j].B_y+u[i][j].B_z*u[i][j].B_z);

            double alpha=0;
            double a_input=1.0,b_input=0,c_input=-1.0;
            u[i][j].rho = (a_input*x+b_input*y+c_input<0?1.0:0.125);
            u[i][j].momentum_x = 0;
            u[i][j].momentum_y = 0;
            u[i][j].momentum_z = 0;
            u[i][j].pressure=(a_input*x+b_input*y+c_input<0?1.0:0.1);
            double shockTube_Bx=0.75,shockTube_By=(a_input*x+b_input*y+c_input<0?1.0:-1.0);
            u[i][j].B_x=shockTube_Bx*cos(alpha)+shockTube_By* cos(alpha+0.5*M_PI);
            u[i][j].B_y=shockTube_Bx* sin(alpha)+shockTube_By* sin(alpha+0.5*M_PI);
            u[i][j].B_z=0;
            u[i][j].energy =u[i][j].pressure/(C_gamma-1)
                             +0.5*(u[i][j].momentum_x*u[i][j].momentum_x+u[i][j].momentum_y*u[i][j].momentum_y+u[i][j].momentum_z*u[i][j].momentum_z)/u[i][j].rho
                             +0.5*(u[i][j].B_x*u[i][j].B_x+u[i][j].B_y*u[i][j].B_y+u[i][j].B_z*u[i][j].B_z);

//            u[i][j].rho = ((x+y)<2?1.0:0.125);
//            u[i][j].momentum_x = 0;
//            u[i][j].momentum_y = 0;
//            u[i][j].momentum_z = 0;
//            u[i][j].energy = ((x+y)<2?1.0:0.1)+0.78125;
//            u[i][j].B_x=0.53033-0.70710*((x+y)<2?1.0:-1.0);
//            u[i][j].B_y=0.53033+0.70710*((x+y)<2?1.0:-1.0);
//            u[i][j].B_z=0;
        }
    }
    uPlus1=u;
    rigidBodyBC(u);
    transmissive(u,NUM_GROW);
//    diagonal_transmissive(u,NUM_GROW);
    double t = tStart;
// 初始化结束

    while (t < tStop) {
        dt = ComputeTimeStep(u, dx, dy, 0.4,NUM_GROW);
        if (t + dt > tStop) {
            dt = tStop - t;
            t = tStop;
//        } else if (t<0.5&&(t+dt)>=0.5){
//            dt=0.5-t;
//            t=0.5;
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
        u=uPlus1;
        rigidBodyBC(u);
        transmissive(u,NUM_GROW);
//        diagonal_transmissive(u,NUM_GROW);
        //end of y update


        //x update
        getFlux(flux, u, dx, dt, 'x',NUM_GROW);
        for (int i = 0; i < x_nCells; i++) {
            for (int j = 0; j < y_nCells; j++) {
                uPlus1[NUM_GROW+i][NUM_GROW+j] = u[NUM_GROW+i][NUM_GROW+j] - (dt / dx) * (flux[i+1][j] - flux[i][j]);
            }
        }
        u=uPlus1;
        rigidBodyBC(u);
        transmissive(u,NUM_GROW);
//        diagonal_transmissive(u,NUM_GROW);
        //end of x update


        //beginning of y update
        getFlux(flux, u, dy, dt/2,'y',NUM_GROW);

        for (int i = 0; i < x_nCells; i++) {
            for (int j = 0; j < y_nCells; j++) {
                uPlus1[NUM_GROW+i][NUM_GROW+j] = u[NUM_GROW+i][NUM_GROW+j] - (dt/2 / dy) * (flux[i][j+1] - flux[i][j]);
            }
        }
        u=uPlus1;
        rigidBodyBC(u);
        transmissive(u,NUM_GROW);
//        diagonal_transmissive(u,NUM_GROW);
        //end of y update

        //coordinate source update//TODO:coordinate
//        for (int i = 0; i < x_nCells; i++) {
//            for (int j = 0; j < y_nCells; j++) {
//                double x = x0 + (i + 0.5-NUM_GROW) * dx;
//                uPlus1[NUM_GROW+i][NUM_GROW+j] = u[NUM_GROW+i][NUM_GROW+j] + dt* u[NUM_GROW+i][NUM_GROW+j].CoordinateSource(1,x);
//            }
//        }
//        transmissive(uPlus1,NUM_GROW);
//        u=uPlus1;

//        DivergencdCleaning(psi,u,psi_flux,dx,dy,dt,NUM_GROW);//TODO:div

        cout<<"time ="<<t<<endl;
    }

    ofstream outfile("\\\\wsl.localhost\\Ubuntu-18.04\\CFD_file\\rigidBody\\Euler.dat");
    if (!outfile)
        return -1;

    for (int i = NUM_GROW; i < x_nCells + NUM_GROW; i++) {
        for (int j = NUM_GROW; j < y_nCells + NUM_GROW; j++) {
            double x = SimMesh.getx(i);
            double y = SimMesh.gety(j);
//            if (get<0>(Phi[i][j])<0) {
//                Equ t(1.0e+5,1.0e+10,1.0e+10,1.0e+10,1.0e+15,1.0e+5,1.0e+5,1.0e+5);
//                outfile << x << " " << y << " " <<t
//                        <<" "<<get<0>(Phi[i][j])
//                        << endl;
//                continue;
//            }
            outfile << x << " " << y << " " << u[i][j]
            <<" "<<get<0>(Phi[i][j])
            << endl;
        }
        outfile<<endl;
    }
    outfile.close();



    ofstream outfile2("\\\\wsl.localhost\\Ubuntu-18.04\\CFD_file\\rigidBody\\Euler_lineout.dat");
    if (!outfile2)
        return -1;

    int resolution=1000;
    double output_dx=(outputEnd_x-outputBegin_x)/(double)resolution;
    double output_dy=(outputEnd_y-outputBegin_y)/(double)resolution;
    for (int i = 0; i < resolution; i++) {
        double x = outputBegin_x+i*output_dx;
        double y = outputBegin_y+i*output_dy;
        double cell_x=(x-x0)/dx-0.5+NUM_GROW;
        double cell_y=(y-y0)/dy-0.5+NUM_GROW;
        int x_1=static_cast<int>(cell_x),y_1=static_cast<int>(cell_y);
        int x_2=x_1+1,y_2=y_1+1;
        Equ result= bilinearInterpolation(u[x_1][y_1],u[x_2][y_1],u[x_1][y_2],u[x_2][y_2],(cell_x-x_1),(cell_y-y_1));
            outfile2 <<(1.5*i/resolution) << " " << result
                    << endl;
    }
    outfile2.close();


//    int tx=10,ty=20;
//    cout<<SimMesh.getx(tx)<<' '<<SimMesh.gety(ty)<<endl;
//    cout<<u[tx][ty][0]<<endl;
//    cout<<get<0>(Phi[tx][ty])<<" the normal is "<<get<1>(Phi[tx][ty]);
//    cout<<"The "<<28<<","<<68<<' '<<"point ";
//    if (  get<0>(Phi[28][68])<=0  &&  get<0>(Phi[28][68])>=(-max(dx,dy)*(NUM_GROW+1))  )
//        for (int k=0;k<8;k++){
//            double T=0.0;
//            cout<<k<<" var :";
//            for (int l = 0; l < Weight[28][68][k].size(); ++l) {
//                tuple<double,double*> t=Weight[28][68][k][l];
//                T=T+get<0>(t)*(*get<1>(t));
//                cout<<"the "<<l<<" is "<<get<0>(t)<<"of"<<*get<1>(t)<<" ";
//            }
//            cout<<endl;
//        }
//    cout<<endl;
    return 0;
}
