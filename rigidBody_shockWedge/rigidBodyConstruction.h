#ifndef Construction_H
#define Construction_H

#include <initializer_list>
#include <vector>
#include <cmath>
#include <tuple>
#include <iostream>

#include "Normal.h"
#include "mesh.h"

using namespace std;

using Pointvalue=tuple<double,double>;

tuple<double,double,double> PointsToLine(Pointvalue P1,Pointvalue P2){
    auto [x1,y1]=P1;
    auto [x2,y2]=P2;
    double a=y1-y2,b=x2-x1,c=x1*y2-x2*y1;
    if (c>0){a*=-1;b*=-1;c*=-1;}
    cout<<a<<' '<<b<<' '<<c<<' '<<endl;
    return make_tuple(a,b,c);
}
double PointsDistance(Pointvalue P1,Pointvalue P2){
    double a=get<0>(P1)-get<0>(P2);
    double b=get<1>(P1)-get<1>(P2);
    return sqrt(a*a+b*b);
}

// 定义基类 Component，作为 func 和 construction 的共同接口
class Component {
public:
    virtual std::tuple<double,Normal> operate(double x, double y) = 0;
    //When making a construction, a principle should be keep. operate return tuple<Phi,Normal_vector>.
    // Phi should always be true and should be minus if it is inside the rigid body.
    // Normal_vector is not always right. it may be adjusted to have correct situation. The right Normal_vector should be unit length point to the boundary no matter outside or inside. (not the gradient of Phi function)
    virtual ~Component() {}
};

// func 的实现，作为 Component 的派生类
class BasicLine : public Component {
    // 使用函数指针而非 std::function
//    double (*func)(double, double);
public:
    double a=0,b=0,c=0;//注意ax+by+c>=0时在模拟区域内
    BasicLine(double at,double bt,double ct):a(at),b(bt),c(ct){}

    Normal getNormal(double x,double y,double Phi){
        double sign = (Phi <= 0) ? 1.0 : -1.0;
        Normal N(sign*a/(sqrt(a*a+b*b)),sign*b/(sqrt(a*a+b*b)));
        return N;
    }

    std::tuple<double,Normal> operate(double x, double y) override {
        double Phi=(a*x+b*y+c)/(sqrt(a*a+b*b));
        Normal N= getNormal(x,y,Phi);
        return std::make_tuple( Phi , N);
    }
};

bool withinLineSegment(BasicLine* e,tuple<double,double> P1,tuple<double,double> P2,double xt,double yt){
    double a=e->a,b=e->b,c=e->c;
    double a1=b,b1=-a,a2=b,b2=-a;
    double x1=get<0>(P1),y1=get<1>(P1);
    double x2=get<0>(P2),y2=get<1>(P2);
    double c1=a*y1-b*x1;
    double c2=a*y2-b*x2;
    if ((a1*x2+b1*y2+c1)<0){a1*=-1;b1*=-1;c1*=-1;}
    if ((a2*x1+b2*y1+c2)<0){a2*=-1;b2*=-1;c2*=-1;}
    bool flag;
//    if ((a1*xt+b1*yt+c1)>=max(SimMesh.dx,SimMesh.dy) && (a2*xt+b2*yt+c2)>=max(SimMesh.dx,SimMesh.dy) && (a*xt+b*yt+c)>=0) flag=true;
    if ((a1*xt+b1*yt+c1)>=0 && (a2*xt+b2*yt+c2)>=0 && (a*xt+b*yt+c)>=0) flag=true;////利用强大的弥合机制 in rigidBodyBC，直接设置为0就完事了
    else flag=false;
    return flag;
}

class BasicCircle : public Component {
public:
    double x=0,y=0,r=0;
    bool inside=true;//模拟区域是inside还是outside
    BasicCircle(double xt,double yt,double rt,bool insideOrFalse):x(xt),y(yt),r(rt),inside(insideOrFalse){}
    BasicCircle(Pointvalue P,double rt,bool insideOrFalse):r(rt),inside(insideOrFalse){
        x=get<0>(P);
        y=get<1>(P);
    }

    Normal getNormal(double a,double b,double Phi){
        double sign = (Phi <= 0) ? 1.0 : -1.0;
        if (Phi!=-r) return Normal(sign*a/(sqrt(a*a+b*b)),sign*b/(sqrt(a*a+b*b))); else return Normal(0,0);
    }

    std::tuple<double,Normal> operate(double xt, double yt) override {
        double a=xt-x;
        double b=yt-y;
        double Phi=sqrt(a*a+b*b)-r;
        Normal N=getNormal(a,b,Phi);
        double sign=(inside)?-1.0:1.0;
        return std::make_tuple( sign*Phi , N);
    }
};

class BasicTriangle : public Component{
public:
    BasicLine *e[3];
    tuple<double,double> coor[3];
    bool inside;//模拟区域是inside还是outside
    BasicTriangle(double x0t,double y0t,double x1t,double y1t,double x2t,double y2t,bool insideOrFalse):inside(insideOrFalse){
        coor[0]= make_tuple(x0t,y0t);
        coor[1]= make_tuple(x1t,y1t);
        coor[2]= make_tuple(x2t,y2t);

        for (int i=0;i<3;i++) {
            int j=i+1; j=j%3;
            auto [a, b, c] = PointsToLine(coor[i], coor[j]);
            int k=j+1; k=k%3;
            if ((a * get<0>(coor[k]) + b * get<1>(coor[k]) + c) > 0) {
                a *= -1;
                b *= -1;
                c *= -1;
            }
            e[i]=new BasicLine(a,b,c);//edge i is made by point i and i+1
        }
    }

    std::tuple<double,Normal> operate(double xt, double yt) override {
        double dist[3];
        Normal N[3];
        for (int i = 0; i <3; ++i) {
            auto [d,n]=e[i]->operate(xt,yt);
            dist[i]=d;
            N[i]=n;
        }

        double T=-1e10;
        Normal NT;
        for (int i=0;i<3;i++){
            if (T<dist[i]) {
                T=dist[i];
                NT = N[i];
            }
        }

        if (dist[0]<0&&dist[1]<0&&dist[2]<0) {////内部点注意离角落很近用center reflected
            Pointvalue currentPoint= make_tuple(xt,yt);
            for (int i=0;i<3;++i){
                double Tdist= PointsDistance(currentPoint,coor[i]);
                if (Tdist<=max(SimMesh.dx,SimMesh.dy)*(SimMesh.NUM_GROW)){
                    BasicCircle cir(coor[i],0,true);
                    auto result=cir.operate(xt,yt);
                    T=get<0>(result);
                    NT=get<1>(result);
                    break;
                }
            }
        }
        else if (withinLineSegment(e[0],coor[0],coor[1],xt,yt)|| withinLineSegment(e[1],coor[1],coor[2],xt,yt)|| withinLineSegment(e[2],coor[2],coor[0],xt,yt)){////正常，算是单边界

        }
        else if ((dist[0]>0&&dist[1]>0)||(dist[1]>0&&dist[2]>0)||(dist[2]>0&&dist[0]>0)){
            for (int i = 0; i < 3; ++i) {
                if (dist[i]<0){
                    int j=(i+1) % 3;
                    int k=(i+2) % 3;
                    BasicCircle cir(coor[k],0,false);
                    auto result=cir.operate(xt,yt);
                    T=get<0>(result);
                    NT=get<1>(result);
                    break;
                }
            }
        }
        else {////三角形的对顶位置和正常位置之间，利用强大的弥合机制就好
            int index=0;
            for (int i=0;i<3;i++){
                if (( (get<0>(coor[index])-xt)*(get<0>(coor[index])-xt)+(get<1>(coor[index])-yt)*(get<1>(coor[index])-yt) ) >  ( (get<0>(coor[i])-xt)*(get<0>(coor[i])-xt)+(get<1>(coor[i])-yt)*(get<1>(coor[i]))-yt) ) {index=i;}
            }
            int index_m1=(index+2) % 3;
            int index_p1=(index+1) % 3;
            Normal N_m1((get<0>(coor[index_m1])-get<0>(coor[index])),(get<1>(coor[index_m1])-get<1>(coor[index])));
            Normal N_p1((get<0>(coor[index_p1])-get<0>(coor[index])),(get<1>(coor[index_p1])-get<1>(coor[index])));
            N_m1=N_m1 *(1/sqrt(N_m1(0)*N_m1(0)*N_m1(1)*N_m1(1))) *(SimMesh.dx/2+SimMesh.dy/2);//TODO:待解决，问题是Phi和normal都不能错，否则。。。暂时这样可能会有错
            N_p1=N_p1 *(1/sqrt(N_p1(0)*N_p1(0)*N_p1(1)*N_p1(1))) *(SimMesh.dx/2+SimMesh.dy/2);
            Normal NNN((get<0>(coor[index])-xt),(get<1>(coor[index])-yt));
            NT=NNN/2+N_m1/4+N_p1/4;
        }


//        if (dist1>0&&dist2>0){
//            BasicCircle cir(x[3],y[3],0,true);
//            auto result=cir.operate(xt,yt);
//            T=get<0>(result);
//            NT=get<1>(result);
//        }
//        if (dist3>0&&dist2>0){
//            BasicCircle cir(x[1],y[1],0,true);
//            auto result=cir.operate(xt,yt);
//            T=get<0>(result);
//            NT=get<1>(result);
//        }
//        if (dist1>0&&dist3>0){
//            BasicCircle cir(x[2],y[2],0,true);
//            auto result=cir.operate(xt,yt);
//            T=get<0>(result);
//            NT=get<1>(result);
//        }

        //before the last step, triangle should always be negative inside.
        double sign=(inside)?-1.0:1.0;
        return make_tuple(sign*T,NT);
    }

    ~BasicTriangle(){
        delete e[0];
        delete e[1];
        delete e[2];
    }
};
//class BasicTriangle : public Component{
//public:
//    BasicLine *e[4];
//    double x[4],y[4];
//    bool inside;
//    BasicTriangle(double x1t,double y1t,double x2t,double y2t,double x3t,double y3t,bool insideOrFalse):inside(insideOrFalse){
//        x[1]=x1t;x[2]=x2t;x[3]=x3t;
//        y[1]=y1t;y[2]=y2t;y[3]=y3t;
//
//        double a3=y[1]-y[2],b3=x[2]-x[1],c3=x[1]*y[2]-x[2]*y[1];
//        if ((a3*x[3]+b3*y[3]+c3)>0){a3*=-1;b3*=-1;c3*=-1;}
//        double a1=y[2]-y[3],b1=x[3]-x[2],c1=x[2]*y[3]-x[3]*y[2];
//        if ((a1*x[1]+b1*y[1]+c1)>0){a1*=-1;b1*=-1;c1*=-1;}
//        double a2=y[3]-y[1],b2=x[1]-x[3],c2=x[3]*y[1]-x[1]*y[3];
//        if ((a2*x[2]+b2*y[2]+c3)>0){a2*=-1;b2*=-1;c2*=-1;}
//
//        e[1]=new BasicLine(a1,b1,c1);
//        e[2]=new BasicLine(a2,b2,c2);
//        e[3]=new BasicLine(a3,b3,c3);
//    }
//
//    std::tuple<double,Normal> operate(double xt, double yt) override {
//        auto [dist1,N1]=e[1]->operate(xt,yt);
//        auto [dist2,N2]=e[2]->operate(xt,yt);
//        auto [dist3,N3]=e[3]->operate(xt,yt);
//
//        double sign=(inside)?-1.0:1.0;
//
//        double T=max({dist1,dist2,dist3});
//        Normal N;if (T==dist1){N=N1;}else if (T==dist2){N=N2;}else {N=N3;}
//
//        if (dist1>0&&dist2>0){
//            BasicCircle cir(x[3],y[3],0,true);
//            auto result=cir.operate(xt,yt);
//            T=get<0>(result);
//            N=get<1>(result);
//        }
//        if (dist3>0&&dist2>0){
//            BasicCircle cir(x[1],y[1],0,true);
//            auto result=cir.operate(xt,yt);
//            T=get<0>(result);
//            N=get<1>(result);
//        }
//        if (dist1>0&&dist3>0){
//            BasicCircle cir(x[2],y[2],0,true);
//            auto result=cir.operate(xt,yt);
//            T=get<0>(result);
//            N=get<1>(result);
//        }
//
//        return make_tuple(sign*T,N);
//    }
//
//    ~BasicTriangle(){
//        delete e[1];
//        delete e[2];
//        delete e[3];
//    }
//};
class BasicRectangle : public Component{
public:
    BasicLine *e[4];
    double coorx[4],coory[4];
    bool inside;
    BasicRectangle(double x1t,double y1t,double x2t,double y2t,double x3t,double y3t,double x4t,double y4t,bool insideOrFalse):inside(insideOrFalse){
        coorx[0]=(x1t);coory[0]=(y1t);
        coorx[1]=(x2t);coory[1]=(y2t);
        coorx[2]=(x3t);coory[2]=(y3t);
        coorx[3]=(x4t);coory[3]=(y4t);

        for (int i=0;i<4;i++){
            int j=i+1;if (j==4) j=0;

            double a=coory[i]-coory[j],b=coorx[j]-coorx[i],c=coorx[i]*coory[j]-coorx[j]*coory[i];

            int it=i-1;if (it==-1) it=3;
            int jt=j+1;if (jt==4) jt=0;

            double dist_it=a*coorx[it]+b*coory[it]+c;
            double dist_jt=a*coorx[jt]+b*coory[jt]+c;

            if (dist_it*dist_jt<0){
                cout<<"Rectangle points sequence is wrongly set.Or it is not a convex rectangle."<<endl;
                exit(-1);
            }

            if (dist_jt>0||dist_it>0){
                a*=-1;
                b*=-1;
                c*=-1;
            }

            e[i]=new BasicLine(a,b,c);
        }
    }

    std::tuple<double,Normal> operate(double x, double y) override {
        auto [dist1,N1]=e[0]->operate(x,y);
        auto [dist2,N2]=e[1]->operate(x,y);
        auto [dist3,N3]=e[2]->operate(x,y);
        auto [dist4,N4]=e[3]->operate(x,y);

        double sign=(inside)?-1.0:1.0;

        double T=max({dist1,dist2,dist3,dist4});
        Normal N;if (T==dist1){N=N1;}else if (T==dist2){N=N2;}else if (T==dist3) {N=N3;} else {N=N4;}

        if (dist1>0&&dist2>0){
            BasicCircle cir(coorx[1],coory[1],0,false);
            auto result=cir.operate(x,y);
            T=get<0>(result);
            N=get<1>(result);
        }
        if (dist2>0&&dist3>0){
            BasicCircle cir(coorx[2],coory[2],0,false);
            auto result=cir.operate(x,y);
            T=get<0>(result);
            N=get<1>(result);
        }
        if (dist3>0&&dist4>0){
            BasicCircle cir(coorx[3],coory[3],0,false);
            auto result=cir.operate(x,y);
            T=get<0>(result);
            N=get<1>(result);
        }
        if (dist4>0&&dist1>0){
            BasicCircle cir(coorx[0],coory[0],0,false);
            auto result=cir.operate(x,y);
            T=get<0>(result);
            N=get<1>(result);
        }


        return make_tuple(sign*T,N);
    }

    ~BasicRectangle(){
        delete e[0];delete e[1];delete e[2];delete e[3];
    }
};

class BasicPolygon:public Component{

};

class Construction : public Component {
public:
    bool and_or = true;
    std::vector<Component*> components;

    // 使用 Component 的指针来初始化
    Construction(std::initializer_list<Component*> list, bool satifyingTrueForAllOrFalseForSome) : and_or(satifyingTrueForAllOrFalseForSome), components(list) {}

    std::tuple<double,Normal> operate(double x, double y) override {
        double sign=(and_or)?1.0:-1.0;
        double T=sign*(1.0e+20);
        Normal N=N_0;
        for (Component*& comp : components){
            std::tuple<double,Normal> result=comp->operate(x,y);
            if ((sign*(get<0>(result))<(sign*T))){T=get<0>(result);N=get<1>(result);}//在角处没有使用扇形处理
        }
        return make_tuple(T,N);
    }

    ~Construction() {
        for (auto& comp : components) {
            delete comp;
        }
    }
};


#endif//boundariesConstruction_H