#ifndef RIGIDBODY_RBMAGNET_H
#define RIGIDBODY_RBMAGNET_H
/**
 * Updating the resistive wall condition with Crank-Nicolson method.
 */

#include "MHDvec.h"
#include "mesh.h"
#include "rigidBodyBC.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include <tuple>
#include <iostream>
#include <cmath>

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> T;
SpMat M;

using namespace std;

Equ magneticCurl(const vector<vector<Equ>> &U,double i,double j){
    double dx = SimMesh.dx;
    double dy = SimMesh.dy;
    Equ T;
    T.B_x=(U[i][j+1].B_z - U[i][j-1].B_z) / (2*dy);
    T.B_y= -(U[i+1][j].B_z - U[i-1][j].B_z) / (2*dx);
    T.B_z=(U[i+1][j].B_y - U[i-1][j].B_y) / (2*dx) - (U[i][j+1].B_x - U[i][j-1].B_x) / (2*dy);
    return T;
}

inline int getidx(int i,int j,int nx,int ny){
    return (j-1)*ny+i;
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

void buildsolver(vector<vector<Equ>> &U,double dt,double eta){
    double dx = SimMesh.dx;
    double dy = SimMesh.dy;
    int x_nCells=SimMesh.x_nCells;
    int y_nCells=SimMesh.y_nCells;
    int NUM_GROW=SimMesh.NUM_GROW;
    int nx=x_nCells+2*NUM_GROW;
    int ny=y_nCells+2*NUM_GROW;

    std::vector<T> triplet_MA,triplet_MB;
    for (int i=0;i<nx;i++){
        for (int j=0;j<ny;j++){
            if (i<NUM_GROW||i>x_nCells+NUM_GROW||j<NUM_GROW||j>y_nCells+NUM_GROW){
                //ghost region first

            } else if (get<0>(Phi[i][j])>0){
                //Plasma area third situation

            } else {
                //rigid bodies second

            }
        }
    }
    SpMat MA,MB;
    MA.setFromTriplets(triplet_MA.begin(), triplet_MA.end());
    MB.setFromTriplets(triplet_MB.begin(), triplet_MB.end());

    Eigen::SparseLU<SpMat> solver;
    solver.analyzePattern(MA);
    solver.factorize(MA);
    M = solver.solve(MB);
}

void updateRigidBodyMagneticField(vector<vector<Equ>> &U,double dt,double eta){
    double dx = SimMesh.dx;
    double dy = SimMesh.dy;
    int x_nCells=SimMesh.x_nCells;
    int y_nCells=SimMesh.y_nCells;
    int NUM_GROW=SimMesh.NUM_GROW;
    int nx=x_nCells+2*NUM_GROW;
    int ny=y_nCells+2*NUM_GROW;

    std::vector<T> triplet_MA,triplet_MB;
    for (int i=0;i<nx;i++){
        for (int j=0;j<ny;j++){

        }
    }
}

#endif//RIGIDBODY_RBMAGNET_H
