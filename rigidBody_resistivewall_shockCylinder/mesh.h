#ifndef RIGIDBODY_MESH_H
#define RIGIDBODY_MESH_H

#include <tuple>
#include "Normal.h"

class mesh{
public:
    int x_nCells,y_nCells,NUM_GROW;
    double x0,x1,y0,y1,dx,dy;
    void set(double x0t,double x1t,int x_nCellst,double y0t,double y1t, int y_nCellst,int NUM_GROWt=1){
        x_nCells=x_nCellst;
        y_nCells=y_nCellst;
        x0=x0t;
        x1=x1t;
        y0=y0t;
        y1=y1t;
        NUM_GROW=NUM_GROWt;

        dx = (x1 - x0) / x_nCells;
        dy = (y1 - y0) / y_nCells;
    }
    Normal get(int i,int j){
        double x = x0 + (i + 0.5-NUM_GROW) * dx;
        double y = y0 + (j + 0.5-NUM_GROW) * dy;

        return Normal (x,y);
    }
    double getx(int i){
        return (x0 + (i + 0.5-NUM_GROW) * dx);
    }
    double gety(int j){
        return (y0 + (j + 0.5-NUM_GROW) * dy);
    }
};

mesh SimMesh;

#endif//RIGIDBODY_MESH_H
