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

BasicCircle* circle1=new BasicCircle(0.9,0.5,0.2, false);
//BasicCircle* circle1=new BasicCircle(1,1,0,true);
Construction rigidBody({circle1}, true);

//BasicLine* line=new BasicLine(-1,0,-1);

//double alpha=45.0    /180.0*M_PI;
//double ax=1.5/2.0*cos(alpha),ay=1.5/2.0*sin(alpha);//length
//double bx=1.0/2.0*cos(alpha+0.5*M_PI),by=1.0/2.0*sin(alpha+0.5*M_PI);//width
//double x_be=1.0,y_be=1.0;
//auto [a_input,b_input,c_input]= PointsToLine(make_tuple(x_be-bx,y_be-by),make_tuple(x_be+bx,y_be+by));
//double outputBegin_x=x_be-ax,outputBegin_y=y_be-ay,outputEnd_x=x_be+ax,outputEnd_y=y_be+ay;
//BasicRectangle* rectangle1=new BasicRectangle(x_be+ax+bx,y_be+ay+by,x_be+ax-bx,y_be+ay-by,x_be-ax-bx,y_be-ay-by,x_be-ax+bx,y_be-ay+by, true);
//Construction rigidBody({rectangle1},false);


//BasicRectangle* rectangle1=new BasicRectangle(0.2,0.6,0.8,0.6,0.8,0.4,0.2,0.4, true);
//Construction rigidBody({rectangle1},true);



//vector<double> valueDBC={0,0,0,0,0,0.75*cos(alpha),0.75*sin(alpha),0},valueNBC(8);//Dirichlet Boundary Condition & Neumann Boundary Condition
vector<double> valueDBC(8),valueNBC(8);
//注意Neumann边界总是对Normal和Tangential求偏导，哪怕是Vx，Vy，Bx，By也是 （而Dirichlet则是x和y方向的值）

//对于不同边界有不同的值，可以建立一个边界值类。()函数类输入位置xy、种类D/N、什么数据int，返回对应边界值类数据的一个引用。

#endif //RIGIDBODY_H