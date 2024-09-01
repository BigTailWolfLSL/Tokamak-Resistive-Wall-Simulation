#ifndef RIGIDBODY_MHDVEC_H
#define RIGIDBODY_MHDVEC_H

#include <vector>
#include <iostream>

#define C_gamma 2.0

class Equ{
public://For IdealMHD gamma=5/3
    double rho=0,momentum_x=0,momentum_y=0,momentum_z=0,energy=0,B_x=0,B_y=0,B_z=0,n=8,
           pressure=0,pressureT=0,gamma=C_gamma,epsilon=0,epsilon_b=0,vel_x=0,vel_y=0,vel_z=0,
           csSqu=0,caSqu=0;
    double& operator[](int nt){
        if (nt==0) return rho;
        if (nt==1) return momentum_x;
        if (nt==2) return momentum_y;
        if (nt==3) return momentum_z;
        if (nt==4) return energy;
        if (nt==5) return B_x;
        if (nt==6) return B_y;
        if (nt==7) return B_z;
        if (nt>=n) exit(1);
    }

    Equ():n(8),gamma(C_gamma){}
    Equ(const Equ& t){
        rho=t.rho;
        momentum_x=t.momentum_x;
        momentum_y=t.momentum_y;
        momentum_z=t.momentum_z;
        energy=t.energy;
        B_x=t.B_x;
        B_y=t.B_y;
        B_z=t.B_z;
        pressure=t.pressure;
        pressureT=t.pressureT;
        gamma=t.gamma;
        epsilon=t.epsilon;
        epsilon_b=t.epsilon_b;
        vel_y=t.vel_y;
        vel_x=t.vel_x;
        vel_z=t.vel_z;
        csSqu=t.csSqu;
        caSqu=t.caSqu;
        n=t.n;
    }
    Equ(double rt,double m_xt,double m_yt,double m_zt,double et,double bxt,double byt,double bzt): rho(rt), momentum_x(m_xt), momentum_y(m_yt),momentum_z(m_zt),energy(et),B_x(bxt),B_y(byt),B_z(bzt),n(8),gamma(C_gamma){}
    void EoS(){
        epsilon_b=0.5*B_z*B_z+0.5*B_x*B_x+0.5*B_y*B_y;
        epsilon= (energy-0.5*momentum_x*momentum_x/rho-0.5*momentum_y*momentum_y/rho-0.5*momentum_z*momentum_z/rho-epsilon_b)/rho;
        pressure=(gamma-1) * rho * epsilon;//ideal MHD
        pressureT=pressure+0.5*B_z*B_z+0.5*B_x*B_x+0.5*B_y*B_y;
        vel_x=momentum_x/rho;
        vel_y=momentum_y/rho;
        vel_z=momentum_z/rho;
        csSqu=gamma*pressure/rho;
        caSqu=(B_x*B_x+B_y*B_y+B_z*B_z)/rho;
        return;
    }
    Equ x_fluxFunc() {
        EoS();
        return Equ((momentum_x),
                   (vel_x*momentum_x) + pressure+0.5*B_y*B_y+0.5*B_z*B_z-0.5*B_x*B_x,
                   (momentum_y*vel_x-B_x*B_y),
                   ((momentum_z*vel_x)-B_x*B_z),
                   (energy + pressure+0.5*B_x*B_x+0.5*B_y*B_y+0.5*B_z*B_z) * vel_x-B_x*(vel_x*B_x+vel_y*B_y+vel_z*B_z),//TODO
                   0,
                   B_y*vel_x-B_x*vel_y,
                   B_z*vel_x-B_x*vel_z
        );
    }

    Equ y_fluxFunc() {
        EoS();
        return Equ((momentum_y),
                   (momentum_x*vel_y-B_x*B_y),
                   momentum_y * vel_y + pressure-0.5*B_y*B_y+0.5*B_z*B_z+0.5*B_x*B_x,
                   momentum_z*vel_y-B_y*B_z,
                   (energy + pressure+0.5*B_x*B_x+0.5*B_y*B_y+0.5*B_z*B_z) *vel_y-(vel_x*B_x+vel_y*B_y+vel_z*B_z)*B_y,
                   B_x*vel_y-B_y*vel_x,
                   0,
                   B_z*vel_y-B_y*vel_z
        );
    }

    Equ CoordinateSource(double alpha,double r){
        EoS();
        Equ t;
        t.rho=rho*vel_x;
        t.momentum_x=rho*vel_x*vel_x;
        //        t.momentum_z=rho*vel_x*vel_z;
        t.energy=vel_x*(energy+pressure);
        t=(-alpha/r)*t;
        return t;
    }

    friend std::ostream &operator<<(std::ostream &os, Equ &equ) {
        equ.EoS();
        os << equ.rho<< ' '
           << equ.vel_x<< ' '
           << equ.vel_y<< ' '
           << equ.vel_z<<' '
           << equ.pressure<<' '
           <<equ.epsilon<<' '
           <<equ.B_x<<' '
           <<equ.B_y<<' '
           <<equ.B_z<<' ';
        return os;
    }
    // TODO:Equ

    Equ operator+(Equ t) {
        Equ Vec_t;
        for (int i = 0; i < n; i++) {
            Vec_t[i] = (*this)[i] + t[i];
        }
        return Equ(Vec_t);
    }

    Equ operator-(Equ t) {
        Equ Vec_t;
        for (int i = 0; i < n; i++) {
            Vec_t[i] = (*this)[i] - t[i];
        }
        return Vec_t;
    }

    Equ operator*(double t) {
        Equ Vec_t;
        for (int i = 0; i < n; i++) {
            Vec_t[i] = (*this)[i] * t;
        }
        return Vec_t;
    }

    Equ operator*(Equ equ) {
        Equ Vec_t;
        for (int i = 0; i < n; i++) {
            Vec_t[i] = (*this)[i] * equ[i];
        }
        return Vec_t;
    }

    friend Equ operator*(double t, Equ equ) {
        Equ Vec_t;
        for (int i = 0; i < equ.n; i++) {
            Vec_t[i] = equ[i] * t;
        }
        return Vec_t;
    }

};

#endif//RIGIDBODY_MHDVEC_H
