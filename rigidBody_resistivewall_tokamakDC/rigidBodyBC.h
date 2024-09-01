#ifndef RIGIDBODY_BC_H
#define RIGIDBODY_BC_H

#include <vector>
#include <tuple>

#include "Eigen/Dense"
#include "Normal.h"
#include "rigidBody.h"
#include "MHDvec.h"
#include "mesh.h"

using namespace std;

//要注意这个函数与传入的u关联起来了

vector<vector<tuple<double,Normal>>> Phi;

using WeightPoint = vector<vector<tuple<double,double*>>>;
vector<vector<WeightPoint>> Weight;

vector<tuple<double,double*>> congregate(vector<tuple<double,double*>> v1,double a1,vector<tuple<double,double*>> v2,double a2){
    vector<tuple<double,double*>> N;
    for (int i=0;i<v1.size();i++){
        N.push_back(make_tuple(get<0>(v1[i])*a1,get<1>(v1[i])));
    }
    for (int i=0;i<v2.size();i++){
        N.push_back(make_tuple(get<0>(v2[i])*a2,get<1>(v2[i])));
    }
    return N;
}
tuple<vector<tuple<double,double*>>,vector<tuple<double,double*>>> seperate(vector<tuple<double,double*>> v,double a1,double a2){
    vector<tuple<double,double*>> N1;
    vector<tuple<double,double*>> N2;
    for (int i=0;i<v.size();i++){
        N1.push_back(make_tuple(get<0>(v[i])*a1,get<1>(v[i])));
        N2.push_back(make_tuple(get<0>(v[i])*a2,get<1>(v[i])));
    }
    return make_tuple(N1,N2);
}
vector<tuple<double,double*>> multiplier(vector<tuple<double,double*>> v,double a){
    vector<tuple<double,double*>> N;
    for (int i=0;i<v.size();i++){
        N.push_back(make_tuple(get<0>(v[i])*a,get<1>(v[i])));
    }
    return N;
}
vector<tuple<double,double*>> addition(vector<tuple<double,double*>> v1,vector<tuple<double,double*>> v2){
    vector<tuple<double,double*>> N;
    for (int i = 0; i < v1.size(); ++i) {
        N.push_back(v1[i]);
    }
    for (int i = 0; i < v2.size(); ++i) {
        N.push_back(v2[i]);
    }
    return N;
}

Normal getBoundaryCoor(int i,int j){
    //get the cloest boundary's coordinate
//    return Normal(0,1);
    return SimMesh.get(i,j)-get<0>(Phi[i][j])*get<1>(Phi[i][j]);//
}

/**
 * 为u[i][j]计算关联更新的边界cells
 * @param U
 * @param i
 * @param j
 */
void getWeight(vector<vector<Equ>> &U,vector<vector<Equ>> &URigidBody,int i,int j){
//The boundary condition for Perfect conducting wall with no ghost fluid way.
    double slipOrNot=1.0;//1 for slip;-1 for no-slip

    if (get<0>(Phi[i][j])<0) {
        Normal cell=Normal(i,j)-Normal(2*get<0>(Phi[i][j])*get<1>(Phi[i][j])(0)/SimMesh.dx,2*get<0>(Phi[i][j])*get<1>(Phi[i][j])(1)/SimMesh.dy);
        Normal coor =SimMesh.get(i,j)-2*get<0>(Phi[i][j])*get<1>(Phi[i][j]);
        Eigen::Vector4d X(1, coor(0),coor(1),coor(0)*coor(1));
        int x1 = static_cast<int>(cell(0));
        int y1 = static_cast<int>(cell(1));
        int x2 = x1 + 1;
        int y2 = y1 + 1;
        Eigen::Matrix<int, 4, 2> bili_cells;
        bili_cells << x1, y1,
                x2, y1,
                x1, y2,
                x2, y2;

        Eigen::Matrix4d A;
        A<< 1 , SimMesh.getx(bili_cells(0,0)) , SimMesh.gety(bili_cells(0,1)) , SimMesh.getx(bili_cells(0,0))*SimMesh.gety(bili_cells(0,1)) ,
                1 , SimMesh.getx(bili_cells(1,0)) , SimMesh.gety(bili_cells(1,1)) , SimMesh.getx(bili_cells(1,0))*SimMesh.gety(bili_cells(1,1)) ,
                1 , SimMesh.getx(bili_cells(2,0)) , SimMesh.gety(bili_cells(2,1)) , SimMesh.getx(bili_cells(2,0))*SimMesh.gety(bili_cells(2,1)) ,
                1 , SimMesh.getx(bili_cells(3,0)) , SimMesh.gety(bili_cells(3,1)) , SimMesh.getx(bili_cells(3,0))*SimMesh.gety(bili_cells(3,1)) ;
        Eigen::Matrix4d A_Dirichlet=A;
        Eigen::Matrix4d A_Neumann=A;
        for (int k=0;k<4;k++){
            if (get<0>(Phi[bili_cells(k,0)][bili_cells(k,1)])<=0){
                Normal IPcoor= getBoundaryCoor(bili_cells(k,0),bili_cells(k,1));
                A_Dirichlet(k,1)=IPcoor(0);
                A_Dirichlet(k,2)=IPcoor(1);
                A_Dirichlet(k,3)=IPcoor(0)*IPcoor(1);
                A_Neumann(k,0)=0;
                double nx=get<1>(Phi[bili_cells(k,0)][bili_cells(k,1)])(0);
                double ny=get<1>(Phi[bili_cells(k,0)][bili_cells(k,1)])(1);
                A_Neumann(k,1)=nx;
                A_Neumann(k,2)=ny;
                A_Neumann(k,3)=nx*IPcoor(1)+ny*IPcoor(0);
            }
        }
        Weight[i][j].resize(U[i][j].n);
        double nx=get<1>(Phi[i][j])(0);//cos theta
        double ny=get<1>(Phi[i][j])(1);//sin theta
        Eigen::Vector4d TD=(A_Dirichlet.inverse().transpose())*X;
        Eigen::Vector4d TN=(A_Neumann.inverse().transpose())*X;

        if (isnan(TD(2))||isnan(TN(2))){//TODO:一个异常强大的弥合机制，在角落处如果interpolation遇到三个或以上的ghost cells变nan时，自动寻找最近的有效点，强制设置为它。
            vector<int> findx,findy;                            //TODO:建议不要随便打开这个if
            int grow=0;
            while (true){
                findx.push_back(x1+grow);findx.push_back(x2+grow);
                findy.push_back(y1+grow);findy.push_back(y2+grow);
                grow++;

                for (int fi=0;fi<findx.size();fi++){
                    for (int fj=0;fj<findy.size();fj++){
                        if (get<0>(Phi[findx[fi]][findy[fj]])>0){
                            for (int fk=0;fk<U[i][j].n;fk++){
                                Weight[i][j][fk].push_back(make_tuple(1,&U[findx[fi]][findy[fj]][fk]));
                            }
                            return;
                        }
                    }
                }
            }
        }


        for (int k=0;k<U[i][j].n;k++){
            if (k==1||k==5) {
                vector<tuple<double,double*>> UxD;UxD.push_back(make_tuple(1,&URigidBody[i][j][k]));
                vector<tuple<double,double*>> UyD;UyD.push_back(make_tuple(1,&URigidBody[i][j][k+1]));
                vector<tuple<double,double*>> UnD= congregate(UxD,nx,UyD,ny);
                vector<tuple<double,double*>> UxN;UxN.push_back(make_tuple(1,&valueNBC[k]));
                vector<tuple<double,double*>> UyN;UyN.push_back(make_tuple(1,&valueNBC[k+1]));
                vector<tuple<double,double*>> UtN= congregate(UxN,ny,UyN,-nx);

                vector<vector<tuple<double,double*>>> vn(4),vt(4);
                for (int l=0;l<4;l++){
                    if (get<0>(Phi[bili_cells(l,0)][bili_cells(l,1)])<=0){

                        vector<tuple<double,double*>> UxD_local;UxD_local.push_back(make_tuple(1,&URigidBody[bili_cells(l,0)][bili_cells(l,1)][k]));
                        vector<tuple<double,double*>> UyD_local;UyD_local.push_back(make_tuple(1,&URigidBody[bili_cells(l,0)][bili_cells(l,1)][k+1]));
                        vector<tuple<double,double*>> UnD_local= congregate(UxD_local,nx,UyD_local,ny);

                        vn[l]=UnD_local;
                        vt[l]=UtN;
                    }
                    else {
                        vn[l].push_back(make_tuple(nx,&U[bili_cells(l,0)][bili_cells(l,1)][k]));//vn=cos(theta)*vxD
                        vn[l].push_back(make_tuple(ny,&U[bili_cells(l,0)][bili_cells(l,1)][k+1]));//+sin(theta)*vyD;

                        vt[l].push_back(make_tuple(ny,&U[bili_cells(l,0)][bili_cells(l,1)][k]));//vt=sin(theta)*vxN
                        vt[l].push_back(make_tuple(-nx,&U[bili_cells(l,0)][bili_cells(l,1)][k+1]));//-cos(theta)vyN;
                    }
                }

                vector<tuple<double,double*>> vn_ip,vt_ip;
                for (int l=0;l<4;l++){
                    vn_ip=addition(vn_ip, multiplier(vn[l],TD(l)));
                    vt_ip=addition(vt_ip, multiplier(vt[l],TN(l)));
                }

                vector<tuple<double,double*>> vn_p=multiplier(vn_ip,-1);
                vector<tuple<double,double*>> vt_p=multiplier(vt_ip,slipOrNot);//TODO:在noslip的时候，应该要考虑Ut
                tuple<vector<tuple<double,double*>>,vector<tuple<double,double*>>> vn_p_pair= seperate(vn_p,nx,ny);
                tuple<vector<tuple<double,double*>>,vector<tuple<double,double*>>> vt_p_pair= seperate(vt_p,ny,-nx);
                tuple<vector<tuple<double,double*>>,vector<tuple<double,double*>>> Un_pair= seperate(UnD,2*nx,2*ny);

                Weight[i][j][k]= addition(get<0>(vn_p_pair),get<0>(vt_p_pair));
                Weight[i][j][k]= addition(Weight[i][j][k],get<0>(Un_pair));
                Weight[i][j][k+1]= addition(get<1>(vn_p_pair),get<1>(vt_p_pair));
                Weight[i][j][k+1]= addition(Weight[i][j][k+1],get<1>(Un_pair));
            }
            else if (k==2||k==6){
                continue;
            }
            else {
                for (int l=0;l<4;l++){
                    if (get<0>(Phi[bili_cells(l,0)][bili_cells(l,1)])<=0){
                        Weight[i][j][k].push_back(make_tuple(TN(l),&valueNBC[k]));
                    } else {
                        Weight[i][j][k].push_back(make_tuple(TN(l),&U[bili_cells(l,0)][bili_cells(l,1)][k]));
                    }
                }
            }
        }
    } else {
        //Exactly on the boundary
        Weight[i][j].resize(U[i][j].n);
        for (int k=0;k<U[i][j].n;k++){
            if (k==1||k==5){
                double nx=get<1>(Phi[i][j])(0);//cos theta
                double ny=get<1>(Phi[i][j])(1);//sin theta

                vector<tuple<double,double*>> Ux;Ux.push_back(make_tuple(1,&URigidBody[i][j][k]));
                vector<tuple<double,double*>> Uy;Uy.push_back(make_tuple(1,&URigidBody[i][j][k+1]));
                vector<tuple<double,double*>> Un= congregate(Ux,nx,Uy,ny);
                vector<tuple<double,double*>> vx;vx.push_back(make_tuple(1,&U[i][j][k]));
                vector<tuple<double,double*>> vy;vy.push_back(make_tuple(1,&U[i][j][k+1]));
                vector<tuple<double,double*>> vt= congregate(vx,ny,vy,-nx);
                tuple<vector<tuple<double,double*>>,vector<tuple<double,double*>>> Un_pair= seperate(Un,nx,ny);
                tuple<vector<tuple<double,double*>>,vector<tuple<double,double*>>> vt_pair= seperate(vt,ny,-nx);

                Weight[i][j][k]= addition(get<0>(Un_pair), multiplier(get<0>(vt_pair),(slipOrNot+1.0)/2.0));//TODO:应该定义Ut并且根据它来得到，noslip是2Ut-vt，slip是vt
                Weight[i][j][k+1]= addition(get<1>(Un_pair), multiplier(get<1>(vt_pair),(slipOrNot+1.0)/2.0));

            } else if (k==2||k==6) {
                continue;
            } else {
                Weight[i][j][k].push_back(make_tuple(1,&U[i][j][k]));
            }
        }
    }
}
//void getWeight(vector<vector<Equ>> &U,vector<vector<Equ>> &URigidBody,int i,int j){
////The boundary condition for resistive wall the first version with Neumann B_tangential and B_z and Dirichlet B_normal.
//    double slipOrNot=1.0;//1 for slip;-1 for no-slip
//
//    if (get<0>(Phi[i][j])<0) {
//        Normal cell=Normal(i,j)-Normal(2*get<0>(Phi[i][j])*get<1>(Phi[i][j])(0)/SimMesh.dx,2*get<0>(Phi[i][j])*get<1>(Phi[i][j])(1)/SimMesh.dy);
//        Normal coor =SimMesh.get(i,j)-2*get<0>(Phi[i][j])*get<1>(Phi[i][j]);
//        Eigen::Vector4d X(1, coor(0),coor(1),coor(0)*coor(1));
//        int x1 = static_cast<int>(cell(0));
//        int y1 = static_cast<int>(cell(1));
//        int x2 = x1 + 1;
//        int y2 = y1 + 1;
//        Eigen::Matrix<int, 4, 2> bili_cells;
//        bili_cells << x1, y1,
//                x2, y1,
//                x1, y2,
//                x2, y2;
//
//        Eigen::Matrix4d A;
//        A<< 1 , SimMesh.getx(bili_cells(0,0)) , SimMesh.gety(bili_cells(0,1)) , SimMesh.getx(bili_cells(0,0))*SimMesh.gety(bili_cells(0,1)) ,
//                1 , SimMesh.getx(bili_cells(1,0)) , SimMesh.gety(bili_cells(1,1)) , SimMesh.getx(bili_cells(1,0))*SimMesh.gety(bili_cells(1,1)) ,
//                1 , SimMesh.getx(bili_cells(2,0)) , SimMesh.gety(bili_cells(2,1)) , SimMesh.getx(bili_cells(2,0))*SimMesh.gety(bili_cells(2,1)) ,
//                1 , SimMesh.getx(bili_cells(3,0)) , SimMesh.gety(bili_cells(3,1)) , SimMesh.getx(bili_cells(3,0))*SimMesh.gety(bili_cells(3,1)) ;
//        Eigen::Matrix4d A_Dirichlet=A;
//        Eigen::Matrix4d A_Neumann=A;
//        for (int k=0;k<4;k++){
//            if (get<0>(Phi[bili_cells(k,0)][bili_cells(k,1)])<=0){
//                Normal IPcoor= getBoundaryCoor(bili_cells(k,0),bili_cells(k,1));
//                A_Dirichlet(k,1)=IPcoor(0);
//                A_Dirichlet(k,2)=IPcoor(1);
//                A_Dirichlet(k,3)=IPcoor(0)*IPcoor(1);
//                A_Neumann(k,0)=0;
//                double nx=get<1>(Phi[bili_cells(k,0)][bili_cells(k,1)])(0);
//                double ny=get<1>(Phi[bili_cells(k,0)][bili_cells(k,1)])(1);
//                A_Neumann(k,1)=nx;
//                A_Neumann(k,2)=ny;
//                A_Neumann(k,3)=nx*IPcoor(1)+ny*IPcoor(0);
//            }
//        }
//        Weight[i][j].resize(U[i][j].n);
//        double nx=get<1>(Phi[i][j])(0);//cos theta
//        double ny=get<1>(Phi[i][j])(1);//sin theta
//        Eigen::Vector4d TD=(A_Dirichlet.inverse().transpose())*X;
//        Eigen::Vector4d TN=(A_Neumann.inverse().transpose())*X;
//
//        if (isnan(TD(2))||isnan(TN(2))){//TODO:一个异常强大的弥合机制，在角落处如果interpolation遇到三个或以上的ghost cells变nan时，自动寻找最近的有效点，强制设置为它。
//            vector<int> findx,findy;                            //TODO:建议不要随便打开这个if
//            int grow=0;
//            while (true){
//                findx.push_back(x1+grow);findx.push_back(x2+grow);
//                findy.push_back(y1+grow);findy.push_back(y2+grow);
//                grow++;
//
//                for (int fi=0;fi<findx.size();fi++){
//                    for (int fj=0;fj<findy.size();fj++){
//                        if (get<0>(Phi[findx[fi]][findy[fj]])>0){
//                            for (int fk=0;fk<U[i][j].n;fk++){
//                                Weight[i][j][fk].push_back(make_tuple(1,&U[findx[fi]][findy[fj]][fk]));
//                            }
//                            return;
//                        }
//                    }
//                }
//            }
//        }
//
//
//        for (int k=0;k<U[i][j].n;k++){
//            if (k==1) {
//                vector<tuple<double,double*>> UxD;UxD.push_back(make_tuple(1,&URigidBody[i][j][k]));
//                vector<tuple<double,double*>> UyD;UyD.push_back(make_tuple(1,&URigidBody[i][j][k+1]));
//                vector<tuple<double,double*>> UnD= congregate(UxD,nx,UyD,ny);
//                vector<tuple<double,double*>> UxN;UxN.push_back(make_tuple(1,&valueNBC[k]));
//                vector<tuple<double,double*>> UyN;UyN.push_back(make_tuple(1,&valueNBC[k+1]));
//                vector<tuple<double,double*>> UtN= congregate(UxN,ny,UyN,-nx);
//
//                vector<vector<tuple<double,double*>>> vn(4),vt(4);
//                for (int l=0;l<4;l++){
//                    if (get<0>(Phi[bili_cells(l,0)][bili_cells(l,1)])<=0){
//                        vector<tuple<double,double*>> UxD_local;UxD_local.push_back(make_tuple(1,&URigidBody[bili_cells(l,0)][bili_cells(l,1)][k]));
//                        vector<tuple<double,double*>> UyD_local;UyD_local.push_back(make_tuple(1,&URigidBody[bili_cells(l,0)][bili_cells(l,1)][k+1]));
//                        vector<tuple<double,double*>> UnD_local= congregate(UxD_local,nx,UyD_local,ny);
//
//                        vn[l]=UnD_local;
//                        vt[l]=UtN;
//                    }
//                    else {
//                        vn[l].push_back(make_tuple(nx,&U[bili_cells(l,0)][bili_cells(l,1)][k]));//vn=cos(theta)*vxD
//                        vn[l].push_back(make_tuple(ny,&U[bili_cells(l,0)][bili_cells(l,1)][k+1]));//+sin(theta)*vyD;
//
//                        vt[l].push_back(make_tuple(ny,&U[bili_cells(l,0)][bili_cells(l,1)][k]));//vt=sin(theta)*vxN
//                        vt[l].push_back(make_tuple(-nx,&U[bili_cells(l,0)][bili_cells(l,1)][k+1]));//-cos(theta)vyN;
//                    }
//                }
//
//                vector<tuple<double,double*>> vn_ip,vt_ip;
//                for (int l=0;l<4;l++){
//                    vn_ip=addition(vn_ip, multiplier(vn[l],TD(l)));
//                    vt_ip=addition(vt_ip, multiplier(vt[l],TN(l)));
//                }
//
//                vector<tuple<double,double*>> vn_p=multiplier(vn_ip,-1);
//                vector<tuple<double,double*>> vt_p=multiplier(vt_ip,slipOrNot);//TODO:在noslip的时候，应该要考虑Ut
//                tuple<vector<tuple<double,double*>>,vector<tuple<double,double*>>> vn_p_pair= seperate(vn_p,nx,ny);
//                tuple<vector<tuple<double,double*>>,vector<tuple<double,double*>>> vt_p_pair= seperate(vt_p,ny,-nx);
//                tuple<vector<tuple<double,double*>>,vector<tuple<double,double*>>> Un_pair= seperate(UnD,2*nx,2*ny);
//
//                Weight[i][j][k]= addition(get<0>(vn_p_pair),get<0>(vt_p_pair));
//                Weight[i][j][k]= addition(Weight[i][j][k],get<0>(Un_pair));
//                Weight[i][j][k+1]= addition(get<1>(vn_p_pair),get<1>(vt_p_pair));
//                Weight[i][j][k+1]= addition(Weight[i][j][k+1],get<1>(Un_pair));
//            }
//            else if (k==2||k==6){
//                continue;
//            }
//            else if (k==5){
//                vector<tuple<double,double*>> UxD;UxD.push_back(make_tuple(1,&URigidBody[i][j][k]));
//                vector<tuple<double,double*>> UyD;UyD.push_back(make_tuple(1,&URigidBody[i][j][k+1]));
//                vector<tuple<double,double*>> UnD= congregate(UxD,nx,UyD,ny);
//                vector<tuple<double,double*>> UxN;UxN.push_back(make_tuple(1,&valueNBC[k]));
//                vector<tuple<double,double*>> UyN;UyN.push_back(make_tuple(1,&valueNBC[k+1]));
//                vector<tuple<double,double*>> UtN= congregate(UxN,ny,UyN,-nx);
//
//                vector<vector<tuple<double,double*>>> vt(4);
//                for (int l=0;l<4;l++){
//                    if (get<0>(Phi[bili_cells(l,0)][bili_cells(l,1)])<=0){
//                        vt[l]=UtN;
//                    }
//                    else {
//                        vt[l].push_back(make_tuple(ny,&U[bili_cells(l,0)][bili_cells(l,1)][k]));//vt=sin(theta)*vxN
//                        vt[l].push_back(make_tuple(-nx,&U[bili_cells(l,0)][bili_cells(l,1)][k+1]));//-cos(theta)vyN;
//                    }
//                }
//
//                vector<tuple<double,double*>> vt_ip;
//                for (int l=0;l<4;l++){
//                    vt_ip=addition(vt_ip, multiplier(vt[l],TN(l)));
//                }
//
//                vector<tuple<double,double*>> vt_p=multiplier(vt_ip,slipOrNot);//TODO:在noslip的时候，应该要考虑Ut
//                tuple<vector<tuple<double,double*>>,vector<tuple<double,double*>>> vt_p_pair= seperate(vt_p,ny,-nx);
//                tuple<vector<tuple<double,double*>>,vector<tuple<double,double*>>> Un_pair= seperate(UnD,nx,ny);
//
//                Weight[i][j][k]= addition(get<0>(Un_pair),get<0>(vt_p_pair));
//                Weight[i][j][k+1]= addition(get<1>(Un_pair),get<1>(vt_p_pair));
//            }
//            else {
//                for (int l=0;l<4;l++){
//                    if (get<0>(Phi[bili_cells(l,0)][bili_cells(l,1)])<=0){
//                        Weight[i][j][k].push_back(make_tuple(TN(l),&valueNBC[k]));
//                    } else {
//                        Weight[i][j][k].push_back(make_tuple(TN(l),&U[bili_cells(l,0)][bili_cells(l,1)][k]));
//                    }
//                }
//            }
//        }
//    } else {
//        //Exactly on the boundary
//        Weight[i][j].resize(U[i][j].n);
//        for (int k=0;k<U[i][j].n;k++){
//            if (k==1||k==5){
//                double nx=get<1>(Phi[i][j])(0);//cos theta
//                double ny=get<1>(Phi[i][j])(1);//sin theta
//
//                vector<tuple<double,double*>> Ux;Ux.push_back(make_tuple(1,&URigidBody[i][j][k]));
//                vector<tuple<double,double*>> Uy;Uy.push_back(make_tuple(1,&URigidBody[i][j][k+1]));
//                vector<tuple<double,double*>> Un= congregate(Ux,nx,Uy,ny);
//                vector<tuple<double,double*>> vx;vx.push_back(make_tuple(1,&U[i][j][k]));
//                vector<tuple<double,double*>> vy;vy.push_back(make_tuple(1,&U[i][j][k+1]));
//                vector<tuple<double,double*>> vt= congregate(vx,ny,vy,-nx);
//                tuple<vector<tuple<double,double*>>,vector<tuple<double,double*>>> Un_pair= seperate(Un,nx,ny);
//                tuple<vector<tuple<double,double*>>,vector<tuple<double,double*>>> vt_pair= seperate(vt,ny,-nx);
//
//                Weight[i][j][k]= addition(get<0>(Un_pair), multiplier(get<0>(vt_pair),(slipOrNot+1.0)/2.0));//TODO:应该定义Ut并且根据它来得到，noslip是2Ut-vt，slip是vt
//                Weight[i][j][k+1]= addition(get<1>(Un_pair), multiplier(get<1>(vt_pair),(slipOrNot+1.0)/2.0));
//
//            } else if (k==2||k==6) {
//                continue;
//            } else {
//                Weight[i][j][k].push_back(make_tuple(1,&U[i][j][k]));
//            }
//        }
//    }
//}
//void getWeight(vector<vector<Equ>> &U,vector<vector<Equ>> &URigidBody,int i,int j){
////The boundary condition for resistive wall rho,momentums,total energy same，B field copy.
//    double slipOrNot=1.0;//1 for slip;-1 for no-slip
//
//    if (get<0>(Phi[i][j])<0) {
//        Normal cell=Normal(i,j)-Normal(2*get<0>(Phi[i][j])*get<1>(Phi[i][j])(0)/SimMesh.dx,2*get<0>(Phi[i][j])*get<1>(Phi[i][j])(1)/SimMesh.dy);
//        Normal coor =SimMesh.get(i,j)-2*get<0>(Phi[i][j])*get<1>(Phi[i][j]);
//        Eigen::Vector4d X(1, coor(0),coor(1),coor(0)*coor(1));
//        int x1 = static_cast<int>(cell(0));
//        int y1 = static_cast<int>(cell(1));
//        int x2 = x1 + 1;
//        int y2 = y1 + 1;
//        Eigen::Matrix<int, 4, 2> bili_cells;
//        bili_cells << x1, y1,
//                x2, y1,
//                x1, y2,
//                x2, y2;
//
//        Eigen::Matrix4d A;
//        A<< 1 , SimMesh.getx(bili_cells(0,0)) , SimMesh.gety(bili_cells(0,1)) , SimMesh.getx(bili_cells(0,0))*SimMesh.gety(bili_cells(0,1)) ,
//                1 , SimMesh.getx(bili_cells(1,0)) , SimMesh.gety(bili_cells(1,1)) , SimMesh.getx(bili_cells(1,0))*SimMesh.gety(bili_cells(1,1)) ,
//                1 , SimMesh.getx(bili_cells(2,0)) , SimMesh.gety(bili_cells(2,1)) , SimMesh.getx(bili_cells(2,0))*SimMesh.gety(bili_cells(2,1)) ,
//                1 , SimMesh.getx(bili_cells(3,0)) , SimMesh.gety(bili_cells(3,1)) , SimMesh.getx(bili_cells(3,0))*SimMesh.gety(bili_cells(3,1)) ;
//        Eigen::Matrix4d A_Dirichlet=A;
//        Eigen::Matrix4d A_Neumann=A;
//        for (int k=0;k<4;k++){
//            if (get<0>(Phi[bili_cells(k,0)][bili_cells(k,1)])<=0){
//                Normal IPcoor= getBoundaryCoor(bili_cells(k,0),bili_cells(k,1));
//                A_Dirichlet(k,1)=IPcoor(0);
//                A_Dirichlet(k,2)=IPcoor(1);
//                A_Dirichlet(k,3)=IPcoor(0)*IPcoor(1);
//                A_Neumann(k,0)=0;
//                double nx=get<1>(Phi[bili_cells(k,0)][bili_cells(k,1)])(0);
//                double ny=get<1>(Phi[bili_cells(k,0)][bili_cells(k,1)])(1);
//                A_Neumann(k,1)=nx;
//                A_Neumann(k,2)=ny;
//                A_Neumann(k,3)=nx*IPcoor(1)+ny*IPcoor(0);
//            }
//        }
//        Weight[i][j].resize(U[i][j].n);
//        double nx=get<1>(Phi[i][j])(0);//cos theta
//        double ny=get<1>(Phi[i][j])(1);//sin theta
//        Eigen::Vector4d TD=(A_Dirichlet.inverse().transpose())*X;
//        Eigen::Vector4d TN=(A_Neumann.inverse().transpose())*X;
//
//        if (isnan(TD(2))||isnan(TN(2))){//TODO:一个异常强大的弥合机制，在角落处如果interpolation遇到三个或以上的ghost cells变nan时，自动寻找最近的有效点，强制设置为它。
//            vector<int> findx,findy;                            //TODO:建议不要随便打开这个if
//            int grow=0;
//            while (true){
//                findx.push_back(x1+grow);findx.push_back(x2+grow);
//                findy.push_back(y1+grow);findy.push_back(y2+grow);
//                grow++;
//
//                for (int fi=0;fi<findx.size();fi++){
//                    for (int fj=0;fj<findy.size();fj++){
//                        if (get<0>(Phi[findx[fi]][findy[fj]])>0){
//                            for (int fk=0;fk<U[i][j].n;fk++){
//                                Weight[i][j][fk].push_back(make_tuple(1,&U[findx[fi]][findy[fj]][fk]));
//                            }
//                            return;
//                        }
//                    }
//                }
//            }
//        }
//
//
//        for (int k=0;k<U[i][j].n;k++){
//            if (k==1) {
//                vector<tuple<double,double*>> UxD;UxD.push_back(make_tuple(1,&URigidBody[i][j][k]));
//                vector<tuple<double,double*>> UyD;UyD.push_back(make_tuple(1,&URigidBody[i][j][k+1]));
//                vector<tuple<double,double*>> UnD= congregate(UxD,nx,UyD,ny);
//                vector<tuple<double,double*>> UxN;UxN.push_back(make_tuple(1,&valueNBC[k]));
//                vector<tuple<double,double*>> UyN;UyN.push_back(make_tuple(1,&valueNBC[k+1]));
//                vector<tuple<double,double*>> UtN= congregate(UxN,ny,UyN,-nx);
//
//                vector<vector<tuple<double,double*>>> vn(4),vt(4);
//                for (int l=0;l<4;l++){
//                    if (get<0>(Phi[bili_cells(l,0)][bili_cells(l,1)])<=0){
//                        vector<tuple<double,double*>> UxD_local;UxD_local.push_back(make_tuple(1,&URigidBody[bili_cells(l,0)][bili_cells(l,1)][k]));
//                        vector<tuple<double,double*>> UyD_local;UyD_local.push_back(make_tuple(1,&URigidBody[bili_cells(l,0)][bili_cells(l,1)][k+1]));
//                        vector<tuple<double,double*>> UnD_local= congregate(UxD_local,nx,UyD_local,ny);
//
//                        vn[l]=UnD_local;
//                        vt[l]=UtN;
//                    }
//                    else {
//                        vn[l].push_back(make_tuple(nx,&U[bili_cells(l,0)][bili_cells(l,1)][k]));//vn=cos(theta)*vxD
//                        vn[l].push_back(make_tuple(ny,&U[bili_cells(l,0)][bili_cells(l,1)][k+1]));//+sin(theta)*vyD;
//
//                        vt[l].push_back(make_tuple(ny,&U[bili_cells(l,0)][bili_cells(l,1)][k]));//vt=sin(theta)*vxN
//                        vt[l].push_back(make_tuple(-nx,&U[bili_cells(l,0)][bili_cells(l,1)][k+1]));//-cos(theta)vyN;
//                    }
//                }
//
//                vector<tuple<double,double*>> vn_ip,vt_ip;
//                for (int l=0;l<4;l++){
//                    vn_ip=addition(vn_ip, multiplier(vn[l],TD(l)));
//                    vt_ip=addition(vt_ip, multiplier(vt[l],TN(l)));
//                }
//
//                vector<tuple<double,double*>> vn_p=multiplier(vn_ip,-1);
//                vector<tuple<double,double*>> vt_p=multiplier(vt_ip,slipOrNot);//TODO:在noslip的时候，应该要考虑Ut
//                tuple<vector<tuple<double,double*>>,vector<tuple<double,double*>>> vn_p_pair= seperate(vn_p,nx,ny);
//                tuple<vector<tuple<double,double*>>,vector<tuple<double,double*>>> vt_p_pair= seperate(vt_p,ny,-nx);
//                tuple<vector<tuple<double,double*>>,vector<tuple<double,double*>>> Un_pair= seperate(UnD,2*nx,2*ny);
//
//                Weight[i][j][k]= addition(get<0>(vn_p_pair),get<0>(vt_p_pair));
//                Weight[i][j][k]= addition(Weight[i][j][k],get<0>(Un_pair));
//                Weight[i][j][k+1]= addition(get<1>(vn_p_pair),get<1>(vt_p_pair));
//                Weight[i][j][k+1]= addition(Weight[i][j][k+1],get<1>(Un_pair));
//            }
//            else if (k==2){
//                continue;
//            }
//            else if (k==5||k==6||k==7){
//                Weight[i][j][k].push_back(make_tuple(1,&URigidBody[i][j][k]));
//            }
//            else {
//                for (int l=0;l<4;l++){
//                    if (get<0>(Phi[bili_cells(l,0)][bili_cells(l,1)])<=0){
//                        Weight[i][j][k].push_back(make_tuple(TN(l),&valueNBC[k]));
//                    } else {
//                        Weight[i][j][k].push_back(make_tuple(TN(l),&U[bili_cells(l,0)][bili_cells(l,1)][k]));
//                    }
//                }
//            }
//        }
//    } else {
//        //Exactly on the boundary
//        Weight[i][j].resize(U[i][j].n);
//        for (int k=0;k<U[i][j].n;k++){
//            if (k==1){
//                double nx=get<1>(Phi[i][j])(0);//cos theta
//                double ny=get<1>(Phi[i][j])(1);//sin theta
//
//                vector<tuple<double,double*>> Ux;Ux.push_back(make_tuple(1,&URigidBody[i][j][k]));
//                vector<tuple<double,double*>> Uy;Uy.push_back(make_tuple(1,&URigidBody[i][j][k+1]));
//                vector<tuple<double,double*>> Un= congregate(Ux,nx,Uy,ny);
//                vector<tuple<double,double*>> vx;vx.push_back(make_tuple(1,&U[i][j][k]));
//                vector<tuple<double,double*>> vy;vy.push_back(make_tuple(1,&U[i][j][k+1]));
//                vector<tuple<double,double*>> vt= congregate(vx,ny,vy,-nx);
//                tuple<vector<tuple<double,double*>>,vector<tuple<double,double*>>> Un_pair= seperate(Un,nx,ny);
//                tuple<vector<tuple<double,double*>>,vector<tuple<double,double*>>> vt_pair= seperate(vt,ny,-nx);
//
//                Weight[i][j][k]= addition(get<0>(Un_pair), multiplier(get<0>(vt_pair),(slipOrNot+1.0)/2.0));//TODO:应该定义Ut并且根据它来得到，noslip是2Ut-vt，slip是vt
//                Weight[i][j][k+1]= addition(get<1>(Un_pair), multiplier(get<1>(vt_pair),(slipOrNot+1.0)/2.0));
//
//            } else if (k==2) {
//                continue;
//            } else if (k==5||k==6||k==7) {
//                Weight[i][j][k].push_back(make_tuple(1,&URigidBody[i][j][k]));
//            } else {
//                Weight[i][j][k].push_back(make_tuple(1,&U[i][j][k]));
//            }
//        }
//    }
//}

/**
 * 带rigidBody的皆是外界调用函数，此处将传入的U和weight绑定起来了
 * @param U is bound into weight. There is a way to just store the i&j in weight but not bound.
 */
void rigidBodyPhi(vector<vector<Equ>> &U,vector<vector<Equ>> &URigidBody){
    double dx = SimMesh.dx;
    double dy = SimMesh.dy;
    int x_nCells=SimMesh.x_nCells;
    int y_nCells=SimMesh.y_nCells;
    int NUM_GROW=SimMesh.NUM_GROW;

    Phi.resize(x_nCells+2*NUM_GROW);
    Weight.resize(x_nCells+2*NUM_GROW);
    for (auto& row : Phi) {
        row.resize(y_nCells+2*NUM_GROW); // 重新定义每行的列数
    }
    for (auto& row : Weight) {
        row.resize(y_nCells+2*NUM_GROW); // 重新定义每行的列数
    }

    for (int i = NUM_GROW; i < x_nCells + NUM_GROW; i++) {
        for (int j = NUM_GROW; j < y_nCells + NUM_GROW; j++) {
            double x = SimMesh.getx(i);
            double y = SimMesh.gety(j);
            Phi[i][j]=rigidBody.operate(x,y);
        }
    }
    for (int i = NUM_GROW; i < x_nCells + NUM_GROW; i++) {
        for (int j = NUM_GROW; j < y_nCells + NUM_GROW; j++) {
            if (  get<0>(Phi[i][j])<=0  &&  get<0>(Phi[i][j])>=(-max(dx,dy)*(NUM_GROW+3))  )  getWeight(U,URigidBody,i,j);
        }
    }
}

/**
 * Update rigid body boundary condition.
 * @param U Update U. should be different to the bound U in rigidBodyPhi.
 */
void rigidBodyBC(vector<vector<Equ>> &U){
    //weight的来源已经被记录下来了，这里的U应当和前面rigidBodyPhi绑定的U不一样
    double dx = SimMesh.dx;
    double dy = SimMesh.dy;
    int x_nCells=SimMesh.x_nCells;
    int y_nCells=SimMesh.y_nCells;
    int NUM_GROW=SimMesh.NUM_GROW;

    for (int i=NUM_GROW;i<x_nCells+NUM_GROW;i++){
        for (int j=NUM_GROW;j<y_nCells+NUM_GROW;j++){
            if (  get<0>(Phi[i][j])<=0  &&  get<0>(Phi[i][j])>=(-max(dx,dy)*(NUM_GROW+3))  )
            for (int k=0;k<U[i][j].n;k++){
                double T=0.0;
                for (int l = 0; l < Weight[i][j][k].size(); ++l) {
                    tuple<double,double*> t=Weight[i][j][k][l];
                    T=T+get<0>(t)*(*get<1>(t));
                }
                U[i][j][k]=T;
            }
        }
    }
}

#endif //RIGIDBODY_BC_H