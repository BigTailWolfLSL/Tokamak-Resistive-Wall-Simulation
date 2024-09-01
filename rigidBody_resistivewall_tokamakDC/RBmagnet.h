#ifndef RIGIDBODY_RBMAGNET_H
#define RIGIDBODY_RBMAGNET_H
//Functions on updating the magnetic field for resistive wall

#include "MHDvec.h"
#include "mesh.h"
#include "rigidBodyBC.h"
#include <Eigen/Dense>
#include <vector>
#include <tuple>
#include <iostream>
#include <boost/numeric/odeint.hpp>
using namespace std;
using namespace boost::numeric::odeint;


Equ magneticCurl(const vector<vector<Equ>> &U,double i,double j){
    double dx = SimMesh.dx;
    double dy = SimMesh.dy;
    Equ T;
    T.B_x=(U[i][j+1].B_z - U[i][j-1].B_z) / (2*dy);
    T.B_y= -(U[i+1][j].B_z - U[i-1][j].B_z) / (2*dx);
    T.B_z=(U[i+1][j].B_y - U[i-1][j].B_y) / (2*dx) - (U[i][j+1].B_x - U[i][j-1].B_x) / (2*dy);
    return T;
}

Equ MCurlTimesCurl(const vector<vector<Equ>> &U,double i,double j){
    double dx = SimMesh.dx;
    double dy = SimMesh.dy;
    Equ T2;
    T2.B_x=1.0/4/dx/dy*(U[i+1][j+1].B_y-U[i-1][j+1].B_y-U[i+1][j-1].B_y+U[i-1][j-1].B_y)
             -1.0/4/dy/dy*(U[i][j+2].B_x-2*U[i][j].B_x+U[i][j-2].B_x);
    T2.B_y=1.0/4/dx/dy*(U[i+1][j+1].B_x-U[i-1][j+1].B_x-U[i+1][j-1].B_x+U[i-1][j-1].B_x)
             -1.0/4/dx/dx*(U[i+2][j].B_y-2*U[i][j].B_y+U[i-2][j].B_y);
    T2.B_z=-1.0/4/dx/dx*(U[i+2][j].B_z-2*U[i][j].B_z+U[i-2][j].B_z)
             -1.0/4/dy/dy*(U[i][j+2].B_z-2*U[i][j].B_z+U[i][j-2].B_z);
    return T2;
}


void setBC(vector<vector<Equ>> &Urigidbody,vector<vector<Equ>> &U){
    double dx = SimMesh.dx;
    double dy = SimMesh.dy;
    int x_nCells=SimMesh.x_nCells;
    int y_nCells=SimMesh.y_nCells;
    int NUM_GROW=SimMesh.NUM_GROW;

    int nx=x_nCells+2*NUM_GROW;
    int ny=y_nCells+2*NUM_GROW;

    for (int i=0;i<nx;i++){
        for (int j=0;j<ny;j++){
            if (get<0>(Phi[i][j])>0){
                Urigidbody[i][j]=U[i][j];
            }
        }
    }
}


void updateRigidBodyMagneticField(vector<vector<Equ>> &U,double dt,double eta){
    double dx = SimMesh.dx;
    double dy = SimMesh.dy;
    int x_nCells=SimMesh.x_nCells;
    int y_nCells=SimMesh.y_nCells;
    int NUM_GROW=SimMesh.NUM_GROW;

    int nx=x_nCells+2*NUM_GROW;
    int ny=y_nCells+2*NUM_GROW;

//    vector<vector<Equ>> Curl(nx,vector<Equ>(ny));
//    for (int i=NUM_GROW;i<NUM_GROW+x_nCells;i++){
//        for (int j=NUM_GROW;j<NUM_GROW+y_nCells;j++){
//            Curl[i][j]= magneticCurl(U,i,j);
//        }
//    }

    vector<vector<Equ>> Uupdate(nx,vector<Equ>(ny));
    for (int i=NUM_GROW;i<NUM_GROW+x_nCells;i++){
        for (int j=NUM_GROW;j<NUM_GROW+y_nCells;j++){
            if (get<0>(Phi[i][j])>0) continue;
            Equ T=MCurlTimesCurl(U,i,j);
            Uupdate[i][j]= U[i][j];
            Uupdate[i][j].B_x=Uupdate[i][j].B_x-dt*eta*T.B_x;
            Uupdate[i][j].B_y=Uupdate[i][j].B_y-dt*eta*T.B_y;
            Uupdate[i][j].B_z=Uupdate[i][j].B_z-dt*eta*T.B_z;
        }
    }
    U=Uupdate;
}

#endif//RIGIDBODY_RBMAGNET_H
