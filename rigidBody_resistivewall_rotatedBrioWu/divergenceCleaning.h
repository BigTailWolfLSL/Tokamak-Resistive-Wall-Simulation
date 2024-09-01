#ifndef RIGIDBODY_DIVERGENCECLEANING_H
#define RIGIDBODY_DIVERGENCECLEANING_H

#include <cmath>
#include "MHDvec.h"
#include "boundary.h"
using namespace std;

/**
 * divergence cleaning
 */
double com_ch(Equ& ul){
    ul.EoS();
    double csSquare=ul.csSqu;
    double caSquare=ul.caSqu;
    double cf_x=sqrt(0.5*(csSquare+caSquare+ sqrt((csSquare+caSquare)*(csSquare+caSquare)-4*csSquare*ul.B_x*ul.B_x/ul.rho)));
    double cf_y=sqrt(0.5*(csSquare+caSquare+ sqrt((csSquare+caSquare)*(csSquare+caSquare)-4*csSquare*ul.B_y*ul.B_y/ul.rho)));
    double cf_z=sqrt(0.5*(csSquare+caSquare+ sqrt((csSquare+caSquare)*(csSquare+caSquare)-4*csSquare*ul.B_z*ul.B_z/ul.rho)));

    double ans;
    ans=max(abs(ul.vel_x)+cf_x, abs(ul.vel_y)+cf_y);
    ans=max(ans,abs(ul.vel_z)+cf_z);
    return ans;
}

void com_divergencd(vector<vector<Equ>> &u, vector<vector<double>> &divergence,double dx,double dy,int NUM_GROW){
    int nx = divergence.size();
    int ny = divergence[0].size();

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            divergence[i][j] = 0;
            divergence[i][j]+=(u[NUM_GROW+i+1][NUM_GROW+j].B_x-u[NUM_GROW+i-1][NUM_GROW+j].B_x)/dx/2;
            divergence[i][j]+=(u[NUM_GROW+i][NUM_GROW+j+1].B_y-u[NUM_GROW+i][NUM_GROW+j-1].B_y)/dy/2;
        }
    }
}

vector<double> DC_xflux(Equ& ul,Equ& ur,double psil,double psir,double ch){
    vector<double> ans(3);
    ans[0]=0.5*ch*(psil-psir)+0.5*ch*ch*(ul.B_x+ur.B_x);
    ans[1]=0.5*ch*(ul.B_x-ur.B_x)+0.5*(psil+psir);
    ans[2]=0;
    return ans;
}

vector<double> DC_yflux(Equ& ul,Equ& ur,double psil,double psir,double ch){
    vector<double> ans(3);
    ans[0]=0.5*ch*(psil-psir)+0.5*ch*ch*(ul.B_y+ur.B_y);
    ans[1]=0;
    ans[2]=0.5*ch*(ul.B_y-ur.B_y)+0.5*(psil+psir);
    return ans;
}

void DivergencdCleaning(vector<vector<double>> &psi, vector<vector<Equ>> &u,vector<vector<vector<double>>> & flux,double dx,double dy,double dt,int NUM_GROW=1){
    int nx = flux.size();
    int ny = flux[0].size();

    double ch=0;
    for (int i=0;i<nx-1;i++){
        for (int j=0;j<ny-1;j++){
            ch=max(ch,com_ch(u[i][j]));
        }
    }

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny-1; j++) {
            flux[i][j]= DC_xflux(u[NUM_GROW+i-1][NUM_GROW+j],u[NUM_GROW+i][NUM_GROW+j],psi[NUM_GROW+i-1][NUM_GROW+j],psi[NUM_GROW+i][NUM_GROW+j],ch);
        }
    }
    for (int i=0;i<nx-1;i++){
        for (int j=0;j<ny-1;j++){
            psi[NUM_GROW+i][NUM_GROW+j] = psi[NUM_GROW+i][NUM_GROW+j] - (dt / dx) * (flux[i+1][j][0] - flux[i][j][0]);
            u[NUM_GROW+i][NUM_GROW+j].B_x = u[NUM_GROW+i][NUM_GROW+j].B_x - (dt / dx) * (flux[i+1][j][1] - flux[i][j][1]);
            u[NUM_GROW+i][NUM_GROW+j].B_y = u[NUM_GROW+i][NUM_GROW+j].B_y - (dt / dx) * (flux[i+1][j][2] - flux[i][j][2]);
        }
    }
    periodic(u,NUM_GROW);
    periodic(psi,NUM_GROW);

    for (int i = 0; i < nx-1; i++) {
        for (int j = 0; j < ny; j++) {
            flux[i][j]= DC_yflux(u[NUM_GROW+i][NUM_GROW+j-1],u[NUM_GROW+i][NUM_GROW+j],psi[NUM_GROW+i][NUM_GROW+j-1],psi[NUM_GROW+i][NUM_GROW+j],ch);
        }
    }
    for (int i=0;i<nx-1;i++){
        for (int j=0;j<ny-1;j++){
            psi[NUM_GROW+i][NUM_GROW+j] = psi[NUM_GROW+i][NUM_GROW+j] - (dt / dy) * (flux[i][j+1][0] - flux[i][j][0]);
            u[NUM_GROW+i][NUM_GROW+j].B_x = u[NUM_GROW+i][NUM_GROW+j].B_x - (dt / dy) * (flux[i][j+1][1] - flux[i][j][1]);
            u[NUM_GROW+i][NUM_GROW+j].B_y = u[NUM_GROW+i][NUM_GROW+j].B_y - (dt / dy) * (flux[i][j+1][2] - flux[i][j][2]);
        }
    }
    periodic(u,NUM_GROW);

    for (int i=0;i<nx-1;i++){
        for (int j=0;j<ny-1;j++){
            psi[NUM_GROW+i][NUM_GROW+j] = psi[NUM_GROW+i][NUM_GROW+j] - dt *ch/0.18*psi[NUM_GROW+i][NUM_GROW+j];
        }
    }
    periodic(psi,NUM_GROW);
}

#endif//RIGIDBODY_DIVERGENCECLEANING_H
