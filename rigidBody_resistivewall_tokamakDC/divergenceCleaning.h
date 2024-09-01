#ifndef RIGIDBODY_DIVERGENCECLEANING_H
#define RIGIDBODY_DIVERGENCECLEANING_H

#include <cmath>
#include "MHDvec.h"
#include "boundary.h"
#include "mesh.h"
#include "rigidBodyBC.h"
using namespace std;

    vector<vector<double>> divergence;//div需要和有效的格子数一致
//    vector<vector<double>> psi(x_nCells+ 2*NUM_GROW,vector<double>(y_nCells+ 2*NUM_GROW));
    vector<vector<vector<double>>> psi_flux;//flux只要比有效格子数大一圈即可

    void initiateDC(){
        int x_nCells=SimMesh.x_nCells;
        int y_nCells=SimMesh.y_nCells;
        divergence.resize(x_nCells,vector<double>(y_nCells));
        psi_flux.resize(x_nCells+1,vector<vector<double>>(y_nCells+1,vector<double>(3)));
    }

/**
 * divergence cleaning
 */
double com_ch(Equ& ul){
    ul.EoS();
    double csSquare=ul.csSqu;
    double caSquare=ul.caSqu;
    double cf_x=sqrt(0.5*(csSquare+caSquare+ sqrt((csSquare+caSquare)*(csSquare+caSquare)-4*csSquare*ul.B_x*ul.B_x/ul.rho)));
    double cf_y=sqrt(0.5*(csSquare+caSquare+ sqrt((csSquare+caSquare)*(csSquare+caSquare)-4*csSquare*ul.B_y*ul.B_y/ul.rho)));
//    double cf_z=sqrt(0.5*(csSquare+caSquare+ sqrt((csSquare+caSquare)*(csSquare+caSquare)-4*csSquare*ul.B_z*ul.B_z/ul.rho)));

    double ans;
    ans=max(abs(ul.vel_x)+cf_x, abs(ul.vel_y)+cf_y);
//    ans=max(ans,abs(ul.vel_z)+cf_z);
    return ans;
}

void com_divergencd(vector<vector<Equ>> &u){
    int nx = divergence.size();
    int ny = divergence[0].size();
    double dx=SimMesh.dx;
    double dy=SimMesh.dy;
    int NUM_GROW = SimMesh.NUM_GROW;

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            divergence[i][j] = 0;
            divergence[i][j]+=(u[NUM_GROW+i+1][NUM_GROW+j].B_x-u[NUM_GROW+i-1][NUM_GROW+j].B_x)/dx/2;
            divergence[i][j]+=(u[NUM_GROW+i][NUM_GROW+j+1].B_y-u[NUM_GROW+i][NUM_GROW+j-1].B_y)/dy/2;
        }
    }
}

vector<double> DC_xflux(Equ& ul,Equ& ur,double ch){
    vector<double> ans(3);
    ans[0]=0.5*ch*(ul.psi-ur.psi)+0.5*ch*ch*(ul.B_x+ur.B_x);
    ans[1]=0.5*ch*(ul.B_x-ur.B_x)+0.5*(ul.psi+ur.psi);
    ans[2]=0;
    return ans;
}

vector<double> DC_yflux(Equ& ul,Equ& ur,double ch){
    vector<double> ans(3);
    ans[0]=0.5*ch*(ul.psi-ur.psi)+0.5*ch*ch*(ul.B_y+ur.B_y);
    ans[1]=0;
    ans[2]=0.5*ch*(ul.B_y-ur.B_y)+0.5*(ul.psi+ur.psi);
    return ans;
}

void DivergencdCleaning(vector<vector<Equ>> &u,double dt){
    int nx = psi_flux.size();
    int ny = psi_flux[0].size();
    int NUM_GROW=SimMesh.NUM_GROW;
    double dx=SimMesh.dx;
    double dy=SimMesh.dy;

    double ch=0;
    for (int i=0;i<nx-1;i++){
        for (int j=0;j<ny-1;j++){
            ch=max(ch,com_ch(u[i][j]));
        }
    }

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny-1; j++) {
            psi_flux[i][j]= DC_xflux(u[NUM_GROW+i-1][NUM_GROW+j],u[NUM_GROW+i][NUM_GROW+j],ch);
        }
    }
    for (int i=0;i<nx-1;i++){
        for (int j=0;j<ny-1;j++){
            double changeE;
            changeE=0.5*u[NUM_GROW+i][NUM_GROW+j].B_x*u[NUM_GROW+i][NUM_GROW+j].B_x
                      +0.5*u[NUM_GROW+i][NUM_GROW+j].B_y*u[NUM_GROW+i][NUM_GROW+j].B_y;
            u[NUM_GROW+i][NUM_GROW+j].psi = u[NUM_GROW+i][NUM_GROW+j].psi - (dt / dx) * (psi_flux[i+1][j][0] - psi_flux[i][j][0]);
            u[NUM_GROW+i][NUM_GROW+j].B_x = u[NUM_GROW+i][NUM_GROW+j].B_x - (dt / dx) * (psi_flux[i+1][j][1] - psi_flux[i][j][1]);
            u[NUM_GROW+i][NUM_GROW+j].B_y = u[NUM_GROW+i][NUM_GROW+j].B_y - (dt / dx) * (psi_flux[i+1][j][2] - psi_flux[i][j][2]);
            changeE=changeE-0.5*u[NUM_GROW+i][NUM_GROW+j].B_x*u[NUM_GROW+i][NUM_GROW+j].B_x
                      -0.5*u[NUM_GROW+i][NUM_GROW+j].B_y*u[NUM_GROW+i][NUM_GROW+j].B_y;
            u[NUM_GROW+i][NUM_GROW+j].energy-=changeE;
        }
    }
    transmissive(u,NUM_GROW);
    rigidBodyBC(u);

    for (int i = 0; i < nx-1; i++) {
        for (int j = 0; j < ny; j++) {
            psi_flux[i][j]= DC_yflux(u[NUM_GROW+i][NUM_GROW+j-1],u[NUM_GROW+i][NUM_GROW+j],ch);
        }
    }
    for (int i=0;i<nx-1;i++){
        for (int j=0;j<ny-1;j++){
            double changeE;
            changeE=0.5*u[NUM_GROW+i][NUM_GROW+j].B_x*u[NUM_GROW+i][NUM_GROW+j].B_x
                      +0.5*u[NUM_GROW+i][NUM_GROW+j].B_y*u[NUM_GROW+i][NUM_GROW+j].B_y;
            u[NUM_GROW+i][NUM_GROW+j].psi = u[NUM_GROW+i][NUM_GROW+j].psi - (dt / dy) * (psi_flux[i][j+1][0] - psi_flux[i][j][0]);
            u[NUM_GROW+i][NUM_GROW+j].B_x = u[NUM_GROW+i][NUM_GROW+j].B_x - (dt / dy) * (psi_flux[i][j+1][1] - psi_flux[i][j][1]);
            u[NUM_GROW+i][NUM_GROW+j].B_y = u[NUM_GROW+i][NUM_GROW+j].B_y - (dt / dy) * (psi_flux[i][j+1][2] - psi_flux[i][j][2]);
            changeE=changeE-0.5*u[NUM_GROW+i][NUM_GROW+j].B_x*u[NUM_GROW+i][NUM_GROW+j].B_x
                      -0.5*u[NUM_GROW+i][NUM_GROW+j].B_y*u[NUM_GROW+i][NUM_GROW+j].B_y;
            u[NUM_GROW+i][NUM_GROW+j].energy-=changeE;
        }
    }
    transmissive(u,NUM_GROW);
    rigidBodyBC(u);

    for (int i=0;i<nx-1;i++){
        for (int j=0;j<ny-1;j++){
            u[NUM_GROW+i][NUM_GROW+j].psi = u[NUM_GROW+i][NUM_GROW+j].psi - dt *ch/0.18*u[NUM_GROW+i][NUM_GROW+j].psi;
        }
    }
    transmissive(u,NUM_GROW);
    rigidBodyBC(u);
}

#endif//RIGIDBODY_DIVERGENCECLEANING_H
