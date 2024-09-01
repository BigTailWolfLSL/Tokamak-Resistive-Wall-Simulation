#include <vector>
#include <tuple>

#include "Normal.h"
#include "rigidBodyConstruction.h"
//#include "MHDvec.h"

using namespace std;

//BasicLine* line1=new BasicLine(-5,0,1);
//BasicLine* line2=new BasicLine(0,-5,1);
//BasicLine* line3=new BasicLine(10,10,-7);
//Construction rigidBody({line1,line2,line3},false);

//BasicTriangle* triangle1=new BasicTriangle(0.2,0.2,0.2,0.5,0.5,0.2, false);
//BasicTriangle* triangle2=new BasicTriangle(0.8,0.8,0.8,0.5,0.5,0.8,false);
//Construction rigidBody({triangle1,triangle2},true);

//BasicTriangle* triangle1=new BasicTriangle(0.2,0.2,0.2,0.5,0.5,0.2, false);
//Construction rigidBody({triangle1},true);

//BasicCircle* circle1=new BasicCircle(0.5,0.5,0.2, false);
//Construction rigidBody({circle1}, false);

//double alpha=30    *(M_PI / 180.0);
//double ax=1.5*cos(alpha),ay=1.5*sin(alpha);//length
//double bx=1.0*cos(alpha+0.5*M_PI),by=1.0*sin(alpha+0.5*M_PI);//width
//double x_be=0.6,y_be=0.2;
//auto [a_input,b_input,c_input]= PointsToLine(make_tuple(x_be+ax/2.0,y_be+ay/2.0),make_tuple(x_be+ax/2.0+bx,y_be+ay/2.0+by));
//double outputBegin_x=x_be+0.5*bx,outputBegin_y=y_be+0.5*by,outputEnd_x=outputBegin_x+ax,outputEnd_y=outputBegin_y+ay;
//BasicRectangle* rectangle1=new BasicRectangle(x_be,y_be,x_be+ax,y_be+ay,x_be+ax+bx,y_be+ay+by,x_be+bx,y_be+by, true);
//Construction rigidBody({rectangle1},false);


//BasicRectangle* rectangle1=new BasicRectangle(0.2,0.6,0.8,0.6,0.8,0.4,0.2,0.4, false);
//Construction rigidBody({rectangle1},true);

BasicTriangle* triangle1=new BasicTriangle(10,0,13,1.555,13,-1.555,false);
Construction rigidBody({triangle1},true);
//BasicTriangle* triangle1=new BasicTriangle(0.10,0,0.13,0.01555,0.13,-0.01555,false);
//Construction rigidBody({triangle1},true);

//BasicCircle* circle1=new BasicCircle(5.5,0,1,false);
//Construction rigidBody({circle1},false);



vector<double> valueDBC(8),valueNBC(8);//Dirichlet Boundary Condition & Neumann Boundary Condition
//注意Neumann边界总是对Normal求偏导，哪怕是Vx，Vy，Bx，By也是
