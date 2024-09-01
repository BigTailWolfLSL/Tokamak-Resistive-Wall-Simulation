#ifndef RIGIDBODY_H
#define RIGIDBODY_H

#include <vector>
#include <tuple>

#include "rigidBodyConstruction.h"
#include "Normal.h"
//#include "MHDvec.h"

using namespace std;

//BasicLine* line1=new BasicLine(-5,0,1);
//BasicLine* line2=new BasicLine(0,-5,1);
//BasicLine* line3=new BasicLine(10,10,-7);
//Construction rigidBody({line1,line2,line3},false);

//BasicTriangle* triangle1=new BasicTriangle(0.2,0.2,0.2,0.5,0.5,0.2, false);
//BasicTriangle* triangle2=new BasicTriangle(0.8,0.8,0.8,0.5,0.5,0.8,false);
//Construction rigidBody({triangle1,triangle2},true);

//BasicTriangle* triangle1=new BasicTriangle(0.2,0.2,0.2,0.5,0.5,0.2, true);
//Construction rigidBody({triangle1},true);

//BasicCircle* circle1=new BasicCircle(0.5,0.5,0.2, false);
//Construction rigidBody({circle1}, false);

//double alpha=45.0    /180.0*M_PI;
//double ax=1.5/2.0*cos(alpha),ay=1.5/2.0*sin(alpha);//length
//double bx=1.0/2.0*cos(alpha+0.5*M_PI),by=1.0/2.0*sin(alpha+0.5*M_PI);//width
//double x_be=1.0,y_be=1.0;
//auto [a_input,b_input,c_input]= PointsToLine(make_tuple(x_be-bx,y_be-by),make_tuple(x_be+bx,y_be+by));
//double outputBegin_x=x_be-ax,outputBegin_y=y_be-ay,outputEnd_x=x_be+ax,outputEnd_y=y_be+ay;
//vector<double> coorx={x_be+ax+bx,x_be+ax-bx,x_be-ax-bx,x_be-ax+bx},coory={y_be+ay+by,y_be+ay-by,y_be-ay-by,y_be-ay+by};


//vector<double> coorx={0.2,0.8,0.8,0.2},coory={0.2,0.2,0.8,0.8};

std::vector<double> coorx = {
         0.5, 0.51, 0.522727, 0.545455, 0.568182, 0.613636, 0.681818, 0.772727, 0.81, 0.863636, 0.93, 1.0, 1.136364, 1.318182, 1.477273, 1.636364, 1.795455, 1.954545, 2.136364, 2.318182, 2.5, 2.659091, 2.772727, 2.863636, 2.931818, 2.977273, 3.0,  3.0, 2.977273, 2.931818, 2.863636, 2.772727, 2.659091, 2.5, 2.318182, 2.136364, 1.954545, 1.795455, 1.636364, 1.477273, 1.318182, 1.136364, 1.0, 0.93, 0.863636, 0.81, 0.772727, 0.681818, 0.613636, 0.568182, 0.545455, 0.522727, 0.51, 0.5
};

std::vector<double> coory = {
         0.973585, 1.086792, 1.2, 1.336000, 1.471698, 1.640189, 1.831321, 2.023, 2.085, 2.150943, 2.21, 2.264151, 2.332075, 2.38, 2.395, 2.377358, 2.332075, 2.263, 2.15, 1.992453, 1.8, 1.59, 1.381132, 1.154717, 0.90566, 0.679245, 0.45283, -0.45283, -0.679245, -0.90566, -1.154717, -1.381132, -1.59, -1.8, -1.992453, -2.15, -2.263, -2.332075, -2.377358, -2.395, -2.38, -2.332075, -2.264151, -2.21, -2.150943, -2.085, -2.023, -1.831321, -1.640189, -1.471698, -1.336000, -1.2, -1.086792, -0.973585
};


BasicConvexPolygon* tokamak=new BasicConvexPolygon(coorx,coory,true);
Construction rigidBody({tokamak},false);



//vector<double> valueDBC={0,0,0,0,0,0.75*cos(alpha),0.75*sin(alpha),0},valueNBC(8);//Dirichlet Boundary Condition & Neumann Boundary Condition
vector<double> valueDBC(8),valueNBC(8);
//注意Neumann边界总是对Normal和Tangential求偏导，哪怕是Vx，Vy，Bx，By也是 （而Dirichlet则是x和y方向的值）

//对于不同边界有不同的值，可以建立一个边界值类。()函数类输入位置xy、种类D/N、什么数据int，返回对应边界值类数据的一个引用。

#endif //RIGIDBODY_H