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

double integratetotalenergy(vector<vector<Equ>> &u){
    int x_ncells=SimMesh.x_nCells;
    int y_ncells=SimMesh.y_nCells;
    double dx=SimMesh.dx;
    double dy=SimMesh.dy;
    int NUM_GROW=SimMesh.NUM_GROW;

    double ans=0;
    for (int i=NUM_GROW;i<x_ncells+NUM_GROW;i++){
        for (int j=NUM_GROW;j<y_ncells+NUM_GROW;j++){
            if (get<0>(Phi[i][j])<=0){ continue;}
            ans+=u[i][j].energy*dx*dy;
        }
    }
    return ans;
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
    double x0 = 0.2;
    double x1 = 3.3;
    double y0 = -2.5;
    double y1 = 2.5;
    int x_nCells = 310;
    int y_nCells = 500;
    double tStart = 0.0;
    double tStop = 0.2;
    double CFL=0.4;
    double checkpointStep=0.002;
    int NUM_GROW=2;//number of ghost cells
    SimMesh.set(x0,x1,x_nCells,y0,y1,y_nCells,NUM_GROW);
    double eta_wall=0.0;

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

            double miu=0.8825;
            double sigma=0.25;

            double xcenter=1.7;
            double ycenter=0;

            double m=2.0;
            double n=1.0;
            double R=10.0;
            double B0=1.0;
            double mu0=4.0*M_PI*1e-7;
            double q0=1.1;

            double J0 = 2.0*B0/(R*mu0*q0);
            double J0_comp = J0 * sqrt(mu0/1e5);
            double constant1=0.055;

            double theta=atan2(y,x);
            double r=sqrt((x-xcenter)*(x-xcenter)+(y-ycenter)*(y-ycenter));
            double Btheta=0.36011444914639353221599919301809*(r<=0.8825?0.2564497898*r:0.1997246991/r);

            u[i][j].rho = (4.95e-1*tanh(20.0*(-r+0.8825))+5.05e-1);
            u[i][j].momentum_x=u[i][j].rho*10*exp(-(r-miu)*(r-miu)/sigma/sigma);
            u[i][j].momentum_y=u[i][j].rho*10*exp(-(r-miu)*(r-miu)/sigma/sigma);
            u[i][j].momentum_z = 0;
            u[i][j].pressure = (4.5e-1*tanh(20.0*(-r+0.8825))+5.5e-1);
            u[i][j].B_x=-Btheta*(y-ycenter)/r;
            u[i][j].B_y=Btheta*(x-xcenter)/r;
            u[i][j].B_z=B0;
            u[i][j].energy =u[i][j].pressure/(C_gamma-1.0)
                                         +0.5*(u[i][j].momentum_x*u[i][j].momentum_x+u[i][j].momentum_y*u[i][j].momentum_y+u[i][j].momentum_z*u[i][j].momentum_z)/u[i][j].rho
                                         +0.5*(u[i][j].B_x*u[i][j].B_x+u[i][j].B_y*u[i][j].B_y+u[i][j].B_z*u[i][j].B_z);

            uRigidBody[i][j].rho=-1;
            uRigidBody[i][j].momentum_x = 0;
            uRigidBody[i][j].momentum_y = 0;
            uRigidBody[i][j].momentum_z = 0;
            uRigidBody[i][j].pressure = 0;
            uRigidBody[i][j].B_x=-Btheta*(y-ycenter)/r;
            uRigidBody[i][j].B_y=Btheta*(x-xcenter)/r;
            uRigidBody[i][j].B_z=B0;
            uRigidBody[i][j].energy =0;
//            u[i][j].energy=(r<=1?1:0);
        }
    }
    uPlus1=u;
    rigidBodyBC(u);
    transmissive(u,NUM_GROW);
    transmissive(uRigidBody,NUM_GROW);

    double t = tStart;
    vector<double> inteTotalEnergy;inteTotalEnergy.push_back(integratetotalenergy(u));
// 初始化结束

    int checkpointi=1;
    while (t < tStop) {
        dt = ComputeTimeStep(u, dx, dy, CFL,NUM_GROW);
        if (t + dt > tStop) {
            dt = tStop - t;
            t = tStop;
        } else if (t<checkpointi*checkpointStep&&(t+dt)>=checkpointStep*checkpointi){
            dt=checkpointStep*checkpointi-t;
            t=checkpointStep*checkpointi;
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

        static const double resistiveSubDt=CFL/eta_wall/2.0/(1/SimMesh.dx/SimMesh.dx+1/SimMesh.dy/SimMesh.dy);
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

        if (t==checkpointStep*checkpointi&&t!=tStop){
            string outfile1="/CFD_file/Tokamak_application/eta00000old/tokamak_t";
            string outfile2=to_string((int)(checkpointi*checkpointStep*1000));
            string outfile3=".dat";
            string outf=outfile1+outfile2+outfile3;
            ofstream outfile(outf);
            if (!outfile)
                return -1;

            inteTotalEnergy.push_back(integratetotalenergy(u));
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
            checkpointi++;
        }
    }
    inteTotalEnergy.push_back(integratetotalenergy(u));
    {
        ofstream outfile("/CFD_file/Tokamak_application/eta00000old/tokamak_t.dat");
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
        ofstream outfile("/CFD_file/Tokamak_application/eta00000old/integrated_totalEnergy.dat");
        if (!outfile)
            return -1;

        for (int i=0;i<inteTotalEnergy.size();i++){
            outfile << i*checkpointStep<<" "
                    << inteTotalEnergy[i]
                    << endl;
        }

        outfile.close();
    }

    return 0;
}
