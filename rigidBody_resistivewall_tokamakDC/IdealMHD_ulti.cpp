#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include <fstream>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <stdio.h>
#include <omp.h>

#include "boundary.h"
#include "MHDvec.h"
#include "mesh.h"
#include "rigidBodyBC.h"
#include "bili_interpolation.h"


using namespace std;

double ComputeTimeStep(vector<vector<Equ>> &u, double dx, double dy,
                       double C_CFL = 1,int NUM_GROW=1) {
    double a = 0;
    int nx = u.size();
    int ny = u[nx - 1].size();
    for (int i = NUM_GROW; i < nx - NUM_GROW; i++) {
        for (int j = NUM_GROW; j < ny - NUM_GROW; j++) {
            u[i][j].EoS();
            double at=max(abs(u[i][j].vel_x)+sqrt(0.5*(u[i][j].csSqu+u[i][j].caSqu+ sqrt((u[i][j].csSqu+u[i][j].caSqu)*(u[i][j].csSqu+u[i][j].caSqu)-4*u[i][j].csSqu*u[i][j].B_x*u[i][j].B_x/u[i][j].rho))),
                          abs(u[i][j].vel_y)+sqrt(0.5*(u[i][j].csSqu+u[i][j].caSqu+ sqrt((u[i][j].csSqu+u[i][j].caSqu)*(u[i][j].csSqu+u[i][j].caSqu)-4*u[i][j].csSqu*u[i][j].B_y*u[i][j].B_y/u[i][j].rho))));
            a=max(a,at);//TODO:DIM 2places to revise
        }
    }
    return C_CFL*min(dx,dy)/a;
}

Equ LaxFriedrichs_flux(Equ &ul, Equ &ur, double const &dx, double const &dt,
                       char direction = 'x') {
// 返回f_i+1/2,ul是i,ur是i+1
    if (direction == 'x') {
        Equ ans = 0.5 * dx / dt * (ul - ur) + 0.5 * (ul.x_fluxFunc() + ur.x_fluxFunc());
        return ans;
    }
    if (direction == 'y') {
        Equ ans = 0.5 * dx / dt * (ul - ur) + 0.5 * (ul.y_fluxFunc() + ur.y_fluxFunc());
        return ans;
    }
    exit(1);
}

Equ Richtmyer_flux(Equ &ul, Equ &ur, double const &dx, double const &dt,
                   char direction = 'x') {
// 返回f_i+1/2,ul是i,ur是i+1
    if (direction == 'x') {
        Equ uMid =
                0.5 * (ul + ur) - 0.5 * dt / dx * (ur.x_fluxFunc() - ul.x_fluxFunc());
        return uMid.x_fluxFunc();
    }
    if (direction == 'y') {
        Equ uMid =
                0.5 * (ul + ur) - 0.5 * dt / dx * (ur.y_fluxFunc() - ul.y_fluxFunc());
        return uMid.y_fluxFunc();
    }
}

Equ FORCE_flux(Equ &ul, Equ &ur, double const &dx, double const &dt,
               char direction = 'x') {
// 返回f_i+1/2,ul是i,ur是i+1
    return 0.5 * (LaxFriedrichs_flux(ul, ur, dx, dt, direction) +
                  Richtmyer_flux(ul, ur, dx, dt, direction));
}

Equ uHLL(Equ &ul, Equ &ur, double const &dx, double const &dt,
         char direction = 'x') {
// // Be careful! it is not flux_HLL but uHLL.
//// For MHD
    ul.EoS();
    ur.EoS();
    if (direction=='x') {
        double csSquare_l=ul.csSqu;
        double csSquare_r=ur.csSqu;
        double caSquare_l=ul.caSqu;
        double caSquare_r=ur.caSqu;
        double cf_l=sqrt(0.5*(csSquare_l+caSquare_l+ sqrt((csSquare_l+caSquare_l)*(csSquare_l+caSquare_l)-4*csSquare_l*ul.B_x*ul.B_x/ul.rho)));
        double cf_r=sqrt(0.5*(csSquare_r+caSquare_r+ sqrt((csSquare_r+caSquare_r)*(csSquare_r+caSquare_r)-4*csSquare_r*ur.B_x*ur.B_x/ur.rho)));
        double sl =min(ul.vel_x,ur.vel_x)-max(cf_l,cf_r);
        double sr =max(ul.vel_x,ur.vel_x)+max(cf_l,cf_r);
        if (sl>0)return ul;
        if (sr<0)return ur;
        if (sl<0&&sr>0){
            Equ uHll=(sr*ur-sl*ul+ul.x_fluxFunc()-ur.x_fluxFunc())*(1/(sr-sl));
            return uHll;
        }
    }
    if (direction=='y'){
        double csSquare_l=ul.csSqu;
        double csSquare_r=ur.csSqu;
        double caSquare_l=ul.caSqu;
        double caSquare_r=ur.caSqu;
        double cf_l=sqrt(0.5*(csSquare_l+caSquare_l+ sqrt((csSquare_l+caSquare_l)*(csSquare_l+caSquare_l)-4*csSquare_l*ul.B_y*ul.B_y/ul.rho)));
        double cf_r=sqrt(0.5*(csSquare_r+caSquare_r+ sqrt((csSquare_r+caSquare_r)*(csSquare_r+caSquare_r)-4*csSquare_r*ur.B_y*ur.B_y/ur.rho)));
        double sl =min(ul.vel_y,ur.vel_y)-max(cf_l,cf_r);
        double sr =max(ul.vel_y,ur.vel_y)+max(cf_l,cf_r);
        if (sl>=0)return ul;
        if (sr<=0)return ur;
        if (sl<0&&sr>0){
            Equ uHll=(sr*ur-sl*ul+ul.x_fluxFunc()-ur.x_fluxFunc())*(1/(sr-sl));
            return uHll;
        }
    }
}

Equ HLL_flux(Equ ul, Equ ur, double const &dx, double const &dt,
             char direction = 'x') {
// // Be careful! it is not flux_HLL but uHLL.
//// For MHD
    ul.EoS();
    ur.EoS();
    if (direction=='x') {
        double csSquare_l=ul.csSqu;
        double csSquare_r=ur.csSqu;
        double caSquare_l=ul.caSqu;
        double caSquare_r=ur.caSqu;
        double cf_l=sqrt(0.5*(csSquare_l+caSquare_l+ sqrt((csSquare_l+caSquare_l)*(csSquare_l+caSquare_l)-4*csSquare_l*ul.B_x*ul.B_x/ul.rho)));
        double cf_r=sqrt(0.5*(csSquare_r+caSquare_r+ sqrt((csSquare_r+caSquare_r)*(csSquare_r+caSquare_r)-4*csSquare_r*ur.B_x*ur.B_x/ur.rho)));
        double sl =min(ul.vel_x,ur.vel_x)-max(cf_l,cf_r);
        double sr =max(ul.vel_x,ur.vel_x)+max(cf_l,cf_r);
        if (sl>=0)return ul.x_fluxFunc();
        if (sr<=0)return ur.x_fluxFunc();
        if (sl<0&&sr>0){
            Equ Flux=(sr*ul.x_fluxFunc()-sl*ur.x_fluxFunc()+sl*sr*(ur-ul)) *(1/(sr-sl));
            return Flux;
        }
    }
    if (direction=='y'){
        double csSquare_l=ul.csSqu;
        double csSquare_r=ur.csSqu;
        double caSquare_l=ul.caSqu;
        double caSquare_r=ur.caSqu;
        double cf_l=sqrt(0.5*(csSquare_l+caSquare_l+ sqrt((csSquare_l+caSquare_l)*(csSquare_l+caSquare_l)-4*csSquare_l*ul.B_y*ul.B_y/ul.rho)));
        double cf_r=sqrt(0.5*(csSquare_r+caSquare_r+ sqrt((csSquare_r+caSquare_r)*(csSquare_r+caSquare_r)-4*csSquare_r*ur.B_y*ur.B_y/ur.rho)));
        double sl =min(ul.vel_y,ur.vel_y)-max(cf_l,cf_r);
        double sr =max(ul.vel_y,ur.vel_y)+max(cf_l,cf_r);
        if (sl>=0)return ul.y_fluxFunc();
        if (sr<=0)return ur.y_fluxFunc();
        if (sl<0&&sr>0){
            Equ Flux=(sr*ul.y_fluxFunc()-sl*ur.y_fluxFunc()+sl*sr*(ur-ul)) *(1/(sr-sl));
            return Flux;
        }
    }
}

Equ HLLC_flux2(Equ &ul, Equ &ur, double const &dx, double const &dt,
               char direction = 'x'){
    ul.EoS();
    ur.EoS();
    if (direction=='x'){
        double csSquare_l=ul.csSqu;
        double csSquare_r=ur.csSqu;
        double caSquare_l=ul.caSqu;
        double caSquare_r=ur.caSqu;
        double cf_l=sqrt(0.5*(csSquare_l+caSquare_l+ sqrt((csSquare_l+caSquare_l)*(csSquare_l+caSquare_l)-4*csSquare_l*ul.B_x*ul.B_x/ul.rho)));
        double cf_r=sqrt(0.5*(csSquare_r+caSquare_r+ sqrt((csSquare_r+caSquare_r)*(csSquare_r+caSquare_r)-4*csSquare_r*ur.B_x*ur.B_x/ur.rho)));
        double sl =min(ul.vel_x,ur.vel_x)-max(cf_l,cf_r);
        double sr =max(ul.vel_x,ur.vel_x)+max(cf_l,cf_r);
        Equ fl=ul.x_fluxFunc();
        Equ fr=ur.x_fluxFunc();//TODO:DIRECTION

        if (sr<0) {
            return ur.x_fluxFunc();
        } else if (sl>0){
            return ul.x_fluxFunc();
        } else {
            Equ ustar=(sr*ur-sl*ul+fl-fr)*(1/(sr-sl));

            double rho_r=ur[0];
            double velx_r=ur[1]/rho_r;
            double vely_r=ur[2]/rho_r;
            double velz_r=ur[3]/rho_r;
            double energy_r=ur[4];
            double Bx_r=ur[5];
            double By_r=ur[6];
            double Bz_r=ur[7];

            double rho_l=ul[0];
            double velx_l=ul[1]/rho_l;
            double vely_l=ul[2]/rho_l;
            double velz_l=ul[3]/rho_l;
            double energy_l=ul[4];
            double Bx_l=ul[5];
            double By_l=ul[6];
            double Bz_l=ul[7];

            double Bx_star=ustar[5];
            double By_star=ustar[6];
            double Bz_star=ustar[7];

            double sstar=ustar[1]/ustar[0];//TODO:DITECTION
            double urbot=velx_r;
            double ulbot=velx_l;
            double brbot=Bx_r;
            double blbot=Bx_l;
            double bstarbot=Bx_star;

            double rhorstar=rho_r*(sr-urbot)/(sr-sstar);
            double rholstar=rho_l*(sl-ulbot)/(sl-sstar);

//            assert(ur.gamma==2);
//            assert(ul.gamma==2);
            double Ptr=(ur.gamma-1)*(energy_r-0.5*rho_r*(velx_r*velx_r+vely_r*vely_r+velz_r*velz_r))+(-ur.gamma+2)*0.5*(ur.B_x*ur.B_x+ur.B_y*ur.B_y+ur.B_z*ur.B_z);
            double Ptl=(ul.gamma-1)*(energy_l-0.5*rho_l*(velx_l*velx_l+vely_l*vely_l+velz_l*velz_l))+(-ul.gamma+2)*0.5*(ul.B_x*ul.B_x+ul.B_y*ul.B_y+ul.B_z*ul.B_z);

            double Ptrstar=Ptr+rho_r*(sr-urbot)*(sstar-urbot)+bstarbot*bstarbot-brbot*brbot;
            double Ptlstar=Ptl+rho_l*(sl-ulbot)*(sstar-ulbot)+bstarbot*bstarbot-blbot*blbot;

            Equ urstar;
            urstar[1]=(ur[1]*(sr-urbot)+Ptrstar-Ptr+brbot*Bx_r-bstarbot*Bx_star)/(sr-sstar);//TODO:DIRECTION
            urstar[2]=(ur[2]*(sr-urbot)            +brbot*By_r-bstarbot*By_star)/(sr-sstar);
            Equ ulstar;
            ulstar[1]=(ul[1]*(sl-ulbot)+Ptlstar-Ptl+blbot*Bx_l-bstarbot*Bx_star)/(sl-sstar);
            ulstar[2]=(ul[2]*(sl-ulbot)            +blbot*By_l-bstarbot*By_star)/(sl-sstar);

            urstar[3]=(ur[3]*(sr-urbot)+brbot*Bz_r-bstarbot*Bz_star)/(sr-sstar);
            ulstar[3]=(ul[3]*(sl-ulbot)+blbot*Bz_l-bstarbot*Bz_star)/(sl-sstar);

            double velx_rstar=urstar[1]/rhorstar;
            double vely_rstar=urstar[2]/rhorstar;
            double velz_rstar=urstar[3]/rhorstar;

            double velx_lstar=ulstar[1]/rholstar;
            double vely_lstar=ulstar[2]/rholstar;
            double velz_lstar=ulstar[3]/rholstar;

            double energy_rstar=((sr-urbot)*energy_r-Ptr*urbot+Ptrstar*sstar+brbot*(velx_r*Bx_r+vely_r*By_r+velz_r*Bz_r)-bstarbot*(velx_rstar*Bx_star+vely_rstar*By_star+velz_rstar*Bz_star))/(sr-sstar);
            double energy_lstar=((sl-ulbot)*energy_l-Ptl*ulbot+Ptlstar*sstar+blbot*(velx_l*Bx_l+vely_l*By_l+velz_l*Bz_l)-bstarbot*(velx_lstar*Bx_star+vely_lstar*By_star+velz_lstar*Bz_star))/(sl-sstar);

            urstar[0]=rhorstar;
            urstar[4]=energy_rstar;
            urstar[5]=ustar[5];
            urstar[6]=ustar[6];
            urstar[7]=ustar[7];

            ulstar[0]=rholstar;
            ulstar[4]=energy_lstar;
            ulstar[5]=ustar[5];
            ulstar[6]=ustar[6];
            ulstar[7]=ustar[7];

            if (sstar>0){
                return fl+sl*(ulstar-ul);
            } else {
                return fr+sr*(urstar-ur);
            }
        }

    } else if (direction=='y'){
        double csSquare_l=ul.csSqu;
        double csSquare_r=ur.csSqu;
        double caSquare_l=ul.caSqu;
        double caSquare_r=ur.caSqu;
        double cf_l=sqrt(0.5*(csSquare_l+caSquare_l+ sqrt((csSquare_l+caSquare_l)*(csSquare_l+caSquare_l)-4*csSquare_l*ul.B_y*ul.B_y/ul.rho)));
        double cf_r=sqrt(0.5*(csSquare_r+caSquare_r+ sqrt((csSquare_r+caSquare_r)*(csSquare_r+caSquare_r)-4*csSquare_r*ur.B_y*ur.B_y/ur.rho)));
        double sl =min(ul.vel_y,ur.vel_y)-max(cf_l,cf_r);
        double sr =max(ul.vel_y,ur.vel_y)+max(cf_l,cf_r);
        Equ fl=ul.y_fluxFunc();
        Equ fr=ur.y_fluxFunc();//TODO:DIRECTION

        if (sr<0) {
            return ur.y_fluxFunc();
        } else if (sl>0){
            return ul.y_fluxFunc();
        } else {
            Equ ustar=(sr*ur-sl*ul+fl-fr)*(1/(sr-sl));

            double rho_r=ur[0];
            double velx_r=ur[1]/rho_r;
            double vely_r=ur[2]/rho_r;
            double velz_r=ur[3]/rho_r;
            double energy_r=ur[4];
            double Bx_r=ur[5];
            double By_r=ur[6];
            double Bz_r=ur[7];

            double rho_l=ul[0];
            double velx_l=ul[1]/rho_l;
            double vely_l=ul[2]/rho_l;
            double velz_l=ul[3]/rho_l;
            double energy_l=ul[4];
            double Bx_l=ul[5];
            double By_l=ul[6];
            double Bz_l=ul[7];

            double Bx_star=ustar[5];
            double By_star=ustar[6];
            double Bz_star=ustar[7];

            double sstar=ustar[2]/ustar[0];//TODO:DITECTION
            double urbot=vely_r;
            double ulbot=vely_l;
            double brbot=By_r;
            double blbot=By_l;
            double bstarbot=By_star;

            double rhorstar=rho_r*(sr-urbot)/(sr-sstar);
            double rholstar=rho_l*(sl-ulbot)/(sl-sstar);

            double Ptr=(ur.gamma-1)*(energy_r-0.5*rho_r*(velx_r*velx_r+vely_r*vely_r+velz_r*velz_r))+(-ur.gamma+2)*0.5*(ur.B_x*ur.B_x+ur.B_y*ur.B_y+ur.B_z*ur.B_z);
            double Ptl=(ul.gamma-1)*(energy_l-0.5*rho_l*(velx_l*velx_l+vely_l*vely_l+velz_l*velz_l))+(-ul.gamma+2)*0.5*(ul.B_x*ul.B_x+ul.B_y*ul.B_y+ul.B_z*ul.B_z);

            double Ptrstar=Ptr+rho_r*(sr-urbot)*(sstar-urbot)+bstarbot*bstarbot-brbot*brbot;
            double Ptlstar=Ptl+rho_l*(sl-ulbot)*(sstar-ulbot)+bstarbot*bstarbot-blbot*blbot;

            Equ urstar;
            urstar[1]=(ur[1]*(sr-urbot)            +brbot*Bx_r-bstarbot*Bx_star)/(sr-sstar);//TODO:DIRECTION
            urstar[2]=(ur[2]*(sr-urbot)+Ptrstar-Ptr+brbot*By_r-bstarbot*By_star)/(sr-sstar);
            Equ ulstar;
            ulstar[1]=(ul[1]*(sl-ulbot)            +blbot*Bx_l-bstarbot*Bx_star)/(sl-sstar);
            ulstar[2]=(ul[2]*(sl-ulbot)+Ptlstar-Ptl+blbot*By_l-bstarbot*By_star)/(sl-sstar);

            urstar[3]=(ur[3]*(sr-urbot)+brbot*Bz_r-bstarbot*Bz_star)/(sr-sstar);
            ulstar[3]=(ul[3]*(sl-ulbot)+blbot*Bz_l-bstarbot*Bz_star)/(sl-sstar);

            double velx_rstar=urstar[1]/rhorstar;
            double vely_rstar=urstar[2]/rhorstar;
            double velz_rstar=urstar[3]/rhorstar;

            double velx_lstar=ulstar[1]/rholstar;
            double vely_lstar=ulstar[2]/rholstar;
            double velz_lstar=ulstar[3]/rholstar;

            double energy_rstar=((sr-urbot)*energy_r-Ptr*urbot+Ptrstar*sstar+brbot*(velx_r*Bx_r+vely_r*By_r+velz_r*Bz_r)-bstarbot*(velx_rstar*Bx_star+vely_rstar*By_star+velz_rstar*Bz_star))/(sr-sstar);
            double energy_lstar=((sl-ulbot)*energy_l-Ptl*ulbot+Ptlstar*sstar+blbot*(velx_l*Bx_l+vely_l*By_l+velz_l*Bz_l)-bstarbot*(velx_lstar*Bx_star+vely_lstar*By_star+velz_lstar*Bz_star))/(sl-sstar);

            urstar[0]=rhorstar;
            urstar[4]=energy_rstar;
            urstar[5]=ustar[5];
            urstar[6]=ustar[6];
            urstar[7]=ustar[7];

            ulstar[0]=rholstar;
            ulstar[4]=energy_lstar;
            ulstar[5]=ustar[5];
            ulstar[6]=ustar[6];
            ulstar[7]=ustar[7];

            if (sstar>0){
                return fl+sl*(ulstar-ul);
            } else {
                return fr+sr*(urstar-ur);
            }
        }
    }
}

Equ uHLLC(Equ &ul, Equ &ur, double const &dx, double const &dt,
          char direction = 'x'){
    ul.EoS();
    ur.EoS();
    if (direction=='x'){
        double csSquare_l=ul.csSqu;
        double csSquare_r=ur.csSqu;
        double caSquare_l=ul.caSqu;
        double caSquare_r=ur.caSqu;
        double cf_l=sqrt(0.5*(csSquare_l+caSquare_l+ sqrt((csSquare_l+caSquare_l)*(csSquare_l+caSquare_l)-4*csSquare_l*ul.B_x*ul.B_x/ul.rho)));
        double cf_r=sqrt(0.5*(csSquare_r+caSquare_r+ sqrt((csSquare_r+caSquare_r)*(csSquare_r+caSquare_r)-4*csSquare_r*ur.B_x*ur.B_x/ur.rho)));
        double sl =min(ul.vel_x,ur.vel_x)-max(cf_l,cf_r);
        double sr =max(ul.vel_x,ur.vel_x)+max(cf_l,cf_r);
        Equ fl=ul.x_fluxFunc();
        Equ fr=ur.x_fluxFunc();//TODO:DIRECTION

        if (sr<0) {
            return ur;
        } else if (sl>0){
            return ul;
        } else {
            Equ ustar=(sr*ur-sl*ul+fl-fr)*(1/(sr-sl));

            double rho_r=ur[0];
            double velx_r=ur[1]/rho_r;
            double vely_r=ur[2]/rho_r;
            double velz_r=ur[3]/rho_r;
            double energy_r=ur[4];
            double Bx_r=ur[5];
            double By_r=ur[6];
            double Bz_r=ur[7];

            double rho_l=ul[0];
            double velx_l=ul[1]/rho_l;
            double vely_l=ul[2]/rho_l;
            double velz_l=ul[3]/rho_l;
            double energy_l=ul[4];
            double Bx_l=ul[5];
            double By_l=ul[6];
            double Bz_l=ul[7];

            double Bx_star=ustar[5];
            double By_star=ustar[6];
            double Bz_star=ustar[7];

            double sstar=ustar[1]/ustar[0];//TODO:DITECTION
            double urbot=velx_r;
            double ulbot=velx_l;
            double brbot=Bx_r;
            double blbot=Bx_l;
            double bstarbot=Bx_star;

            double rhorstar=rho_r*(sr-urbot)/(sr-sstar);
            double rholstar=rho_l*(sl-ulbot)/(sl-sstar);

            assert(ur.gamma==2);
            assert(ul.gamma==2);
            double Ptr=(ur.gamma-1)*(energy_r-0.5*rho_r*(velx_r*velx_r+vely_r*vely_r+velz_r*velz_r))+(-ur.gamma+2)*0.5*(ur.B_x*ur.B_x+ur.B_y*ur.B_y+ur.B_z*ur.B_z);
            double Ptl=(ul.gamma-1)*(energy_l-0.5*rho_l*(velx_l*velx_l+vely_l*vely_l+velz_l*velz_l))+(-ul.gamma+2)*0.5*(ul.B_x*ul.B_x+ul.B_y*ul.B_y+ul.B_z*ul.B_z);

            double Ptrstar=Ptr+rho_r*(sr-urbot)*(sstar-urbot)+bstarbot*bstarbot-brbot*brbot;
            double Ptlstar=Ptl+rho_l*(sl-ulbot)*(sstar-ulbot)+bstarbot*bstarbot-blbot*blbot;

            Equ urstar;
            urstar[1]=(ur[1]*(sr-urbot)+Ptrstar-Ptr+brbot*Bx_r-bstarbot*Bx_star)/(sr-sstar);//TODO:DIRECTION
            urstar[2]=(ur[2]*(sr-urbot)            +brbot*By_r-bstarbot*By_star)/(sr-sstar);
            Equ ulstar;
            ulstar[1]=(ul[1]*(sl-ulbot)+Ptlstar-Ptl+blbot*Bx_l-bstarbot*Bx_star)/(sl-sstar);
            ulstar[2]=(ul[2]*(sl-ulbot)            +blbot*By_l-bstarbot*By_star)/(sl-sstar);

            urstar[3]=(ur[3]*(sr-urbot)+brbot*Bz_r-bstarbot*Bz_star)/(sr-sstar);
            ulstar[3]=(ul[3]*(sl-ulbot)+blbot*Bz_l-bstarbot*Bz_star)/(sl-sstar);

            double velx_rstar=urstar[1]/rhorstar;
            double vely_rstar=urstar[2]/rhorstar;
            double velz_rstar=urstar[3]/rhorstar;

            double velx_lstar=ulstar[1]/rholstar;
            double vely_lstar=ulstar[2]/rholstar;
            double velz_lstar=ulstar[3]/rholstar;

            double energy_rstar=((sr-urbot)*energy_r-Ptr*urbot+Ptrstar*sstar+brbot*(velx_r*Bx_r+vely_r*By_r+velz_r*Bz_r)-bstarbot*(velx_rstar*Bx_star+vely_rstar*By_star+velz_rstar*Bz_star))/(sr-sstar);
            double energy_lstar=((sl-ulbot)*energy_l-Ptl*ulbot+Ptlstar*sstar+blbot*(velx_l*Bx_l+vely_l*By_l+velz_l*Bz_l)-bstarbot*(velx_lstar*Bx_star+vely_lstar*By_star+velz_lstar*Bz_star))/(sl-sstar);

            urstar[0]=rhorstar;
            urstar[4]=energy_rstar;
            urstar[5]=ustar[5];
            urstar[6]=ustar[6];
            urstar[7]=ustar[7];

            ulstar[0]=rholstar;
            ulstar[4]=energy_lstar;
            ulstar[5]=ustar[5];
            ulstar[6]=ustar[6];
            ulstar[7]=ustar[7];

            if (sstar>0){
                return ulstar;
            } else {
                return urstar;
            }
        }

    } else if (direction=='y'){
        double csSquare_l=ul.csSqu;
        double csSquare_r=ur.csSqu;
        double caSquare_l=ul.caSqu;
        double caSquare_r=ur.caSqu;
        double cf_l=sqrt(0.5*(csSquare_l+caSquare_l+ sqrt((csSquare_l+caSquare_l)*(csSquare_l+caSquare_l)-4*csSquare_l*ul.B_y*ul.B_y/ul.rho)));
        double cf_r=sqrt(0.5*(csSquare_r+caSquare_r+ sqrt((csSquare_r+caSquare_r)*(csSquare_r+caSquare_r)-4*csSquare_r*ur.B_y*ur.B_y/ur.rho)));
        double sl =min(ul.vel_y,ur.vel_y)-max(cf_l,cf_r);
        double sr =max(ul.vel_y,ur.vel_y)+max(cf_l,cf_r);
        Equ fl=ul.y_fluxFunc();
        Equ fr=ur.y_fluxFunc();//TODO:DIRECTION

        if (sr<0) {
            return ur;
        } else if (sl>0){
            return ul;
        } else {
            Equ ustar=(sr*ur-sl*ul+fl-fr)*(1/(sr-sl));

            double rho_r=ur[0];
            double velx_r=ur[1]/rho_r;
            double vely_r=ur[2]/rho_r;
            double velz_r=ur[3]/rho_r;
            double energy_r=ur[4];
            double Bx_r=ur[5];
            double By_r=ur[6];
            double Bz_r=ur[7];

            double rho_l=ul[0];
            double velx_l=ul[1]/rho_l;
            double vely_l=ul[2]/rho_l;
            double velz_l=ul[3]/rho_l;
            double energy_l=ul[4];
            double Bx_l=ul[5];
            double By_l=ul[6];
            double Bz_l=ul[7];

            double Bx_star=ustar[5];
            double By_star=ustar[6];
            double Bz_star=ustar[7];

            double sstar=ustar[2]/ustar[0];//TODO:DITECTION
            double urbot=vely_r;
            double ulbot=vely_l;
            double brbot=By_r;
            double blbot=By_l;
            double bstarbot=By_star;

            double rhorstar=rho_r*(sr-urbot)/(sr-sstar);
            double rholstar=rho_l*(sl-ulbot)/(sl-sstar);

            double Ptr=(ur.gamma-1)*(energy_r-0.5*rho_r*(velx_r*velx_r+vely_r*vely_r+velz_r*velz_r))+(-ur.gamma+2)*0.5*(ur.B_x*ur.B_x+ur.B_y*ur.B_y+ur.B_z*ur.B_z);
            double Ptl=(ul.gamma-1)*(energy_l-0.5*rho_l*(velx_l*velx_l+vely_l*vely_l+velz_l*velz_l))+(-ul.gamma+2)*0.5*(ul.B_x*ul.B_x+ul.B_y*ul.B_y+ul.B_z*ul.B_z);

            double Ptrstar=Ptr+rho_r*(sr-urbot)*(sstar-urbot)+bstarbot*bstarbot-brbot*brbot;
            double Ptlstar=Ptl+rho_l*(sl-ulbot)*(sstar-ulbot)+bstarbot*bstarbot-blbot*blbot;

            Equ urstar;
            urstar[1]=(ur[1]*(sr-urbot)            +brbot*Bx_r-bstarbot*Bx_star)/(sr-sstar);//TODO:DIRECTION
            urstar[2]=(ur[2]*(sr-urbot)+Ptrstar-Ptr+brbot*By_r-bstarbot*By_star)/(sr-sstar);
            Equ ulstar;
            ulstar[1]=(ul[1]*(sl-ulbot)            +blbot*Bx_l-bstarbot*Bx_star)/(sl-sstar);
            ulstar[2]=(ul[2]*(sl-ulbot)+Ptlstar-Ptl+blbot*By_l-bstarbot*By_star)/(sl-sstar);

            urstar[3]=(ur[3]*(sr-urbot)+brbot*Bz_r-bstarbot*Bz_star)/(sr-sstar);
            ulstar[3]=(ul[3]*(sl-ulbot)+blbot*Bz_l-bstarbot*Bz_star)/(sl-sstar);

            double velx_rstar=urstar[1]/rhorstar;
            double vely_rstar=urstar[2]/rhorstar;
            double velz_rstar=urstar[3]/rhorstar;

            double velx_lstar=ulstar[1]/rholstar;
            double vely_lstar=ulstar[2]/rholstar;
            double velz_lstar=ulstar[3]/rholstar;

            double energy_rstar=((sr-urbot)*energy_r-Ptr*urbot+Ptrstar*sstar+brbot*(velx_r*Bx_r+vely_r*By_r+velz_r*Bz_r)-bstarbot*(velx_rstar*Bx_star+vely_rstar*By_star+velz_rstar*Bz_star))/(sr-sstar);
            double energy_lstar=((sl-ulbot)*energy_l-Ptl*ulbot+Ptlstar*sstar+blbot*(velx_l*Bx_l+vely_l*By_l+velz_l*Bz_l)-bstarbot*(velx_lstar*Bx_star+vely_lstar*By_star+velz_lstar*Bz_star))/(sl-sstar);

            urstar[0]=rhorstar;
            urstar[4]=energy_rstar;
            urstar[5]=ustar[5];
            urstar[6]=ustar[6];
            urstar[7]=ustar[7];

            ulstar[0]=rholstar;
            ulstar[4]=energy_lstar;
            ulstar[5]=ustar[5];
            ulstar[6]=ustar[6];
            ulstar[7]=ustar[7];

            if (sstar>0){
                return ulstar;
            } else {
                return urstar;
            }
        }
    }
}

double limiter(double r){
    if (r<=0){
        return 0;
    } else {
        return min((2/(1+r)) ,(2*r/(1+r)));//VanLeer
    }
}

Equ MUSCL_Hancock_flux(vector<vector<Equ>> &u,int i,int j, double const &dx, double const &dt,
                       char direction = 'x',double omega=0) {
// 返回f_i+1/2,ul是i,ur是i+1
//// Be really careful that NUM_GROW must be >=2.
//// Using HLLC.
    Equ deltaHead_l;
    Equ deltaHead_r;
    if (direction=='x'){
        {
            Equ delta_l = u[i - 1][j] - u[i - 2][j];
            Equ delta_r = u[i][j] - u[i - 1][j];
            double r;
            for (int l = 0; l < delta_l.n; l++) {
                r = delta_l[l] / (delta_r[l] + 1.0e-20);//(delta_r[l] != 0 ? delta_l[l] / delta_r[l] : 1.0e+20);
                deltaHead_l[l]=0.5*(1+omega)*delta_l[l]+0.5*(1-omega)*delta_r[l];
                deltaHead_l[l]=deltaHead_l[l]* limiter(r);
            }
        }
        {
            Equ delta_l = u[i][j] - u[i - 1][j];
            Equ delta_r = u[i+1][j] - u[i][j];
            double r;
            for (int l = 0; l < delta_l.n; l++) {
                r = delta_l[l] / (delta_r[l] + 1.0e-20);//(delta_r[l] != 0 ? delta_l[l] / delta_r[l] : 1.0e+20);
                deltaHead_r[l]=0.5*(1+omega)*delta_l[l]+0.5*(1-omega)*delta_r[l];
                deltaHead_r[l]=deltaHead_r[l]* limiter(r);//VanLeer
            }
        }

        Equ ul=(u[i-1][j]+0.5*deltaHead_l)-0.5*dt/dx*((u[i-1][j]+0.5*deltaHead_l).x_fluxFunc() -(u[i-1][j]-0.5*deltaHead_l).x_fluxFunc());
        Equ ur=(u[i][j]-0.5*deltaHead_r)-0.5*dt/dx*((u[i][j]+0.5*deltaHead_r).x_fluxFunc() -(u[i][j]-0.5*deltaHead_r).x_fluxFunc());
        return HLLC_flux2(ul,ur,dx,dt,direction);//TODO
    }
    if (direction=='y'){
        {
            Equ delta_l = u[i][j- 1] - u[i][j- 2];
            Equ delta_r = u[i][j] - u[i ][j- 1];
            double r;
            for (int l = 0; l < delta_l.n; l++) {
                r = delta_l[l] / (delta_r[l] + 1.0e-20);//(delta_r[l] != 0 ? delta_l[l] / delta_r[l] : 1.0e+20);
                deltaHead_l[l]=0.5*(1+omega)*delta_l[l]+0.5*(1-omega)*delta_r[l];
                deltaHead_l[l]=deltaHead_l[l]* limiter(r);//VanLeer
            }

        }
        {
            Equ delta_l = u[i][j] - u[i][j-1];
            Equ delta_r = u[i][j+1] - u[i][j];
            double r;
            for (int l = 0; l < delta_l.n; l++) {
                r = delta_l[l] / (delta_r[l] + 1.0e-20);//(delta_r[l] != 0 ? delta_l[l] / delta_r[l] : 1.0e+20);
                deltaHead_r[l]=0.5*(1+omega)*delta_l[l]+0.5*(1-omega)*delta_r[l];
                deltaHead_r[l]=deltaHead_r[l]*limiter(r);//VanLeer
            }
        }

        Equ ul=(u[i][j-1]+0.5*deltaHead_l)-0.5*dt/dx*((u[i][j-1]+0.5*deltaHead_l).y_fluxFunc() -(u[i][j-1]-0.5*deltaHead_l).y_fluxFunc());
        Equ ur=(u[i][j]-0.5*deltaHead_r)-0.5*dt/dx*((u[i][j]+0.5*deltaHead_r).y_fluxFunc() -(u[i][j]-0.5*deltaHead_r).y_fluxFunc());
        return HLLC_flux2(ul,ur,dx,dt,direction);
    }
}

void getFlux(vector<vector<Equ>> &flux, vector<vector<Equ>> &u, double dx,
             double dt, char direction = 'x',int NUM_GROW=1) {
// 'dx' is already the delta in the 'direction'.
    int nx = flux.size();
    int ny = flux[0].size();
    if (direction == 'x') {
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny-1; j++) {
//                flux[i][j] = HLLC_flux2(u[NUM_GROW+i-1][NUM_GROW+j],u[NUM_GROW+i][NUM_GROW+j], dx, dt,direction);
                flux[i][j]= MUSCL_Hancock_flux(u,NUM_GROW+i,NUM_GROW+j,dx, dt,direction);
//TODO:We change the x-direction numerical method here.
            }
        }
    }
    if (direction == 'y') {
        for (int i = 0; i < nx-1; i++) {
            for (int j = 0; j < ny; j++) {
//                flux[i][j] = HLLC_flux2(u[NUM_GROW+i][NUM_GROW+j-1],u[NUM_GROW+i][NUM_GROW+j], dx, dt,direction);
                flux[i][j]= MUSCL_Hancock_flux(u,NUM_GROW+i,NUM_GROW+j,dx, dt,direction);
//TODO:We change the y-direction numerical method here.
            }
        }
    }
    return;
}

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

int main() {
    double x0 = 0.0;
    double x1 = 2.0;
    double y0 = 0.0;
    double y1 = 2.0;
    int x_nCells = 100;
    int y_nCells = 100;
    double tStart = 0.0;
    double tStop = 0.20;
    int NUM_GROW=2;//number of ghost cells
    SimMesh.set(x0,x1,x_nCells,y0,y1,y_nCells,NUM_GROW);

    double dx = SimMesh.dx;
    double dy = SimMesh.dy;
    double dt;
    vector<vector<Equ>> u(x_nCells + 2*NUM_GROW,
                          vector<Equ>(y_nCells + 2*NUM_GROW)); // 有效范围1~nCells+1
    vector<vector<Equ>> uPlus1(x_nCells + 2*NUM_GROW,
                               vector<Equ>(y_nCells + 2*NUM_GROW)); // 两个u之间有一个flux
    rigidBodyPhi(uPlus1);
    vector<vector<Equ>> flux(x_nCells + 1,
                             vector<Equ>(y_nCells + 1)); // flux[i]=f_i+1/2


//    vector<vector<double>> divergence(x_nCells,vector<double>(y_nCells));//TODO:div
//    vector<vector<double>> psi(x_nCells+ 2*NUM_GROW,vector<double>(y_nCells+ 2*NUM_GROW));
//    vector<vector<vector<double>>> psi_flux(x_nCells+1,vector<vector<double>>(y_nCells+1,vector<double>(3)));


// 2D时，flux计算x方向上相当于x方向加了1/2，计算y方向上相当于y方向加了1/2
// 定义结束

    for (int i = NUM_GROW; i < x_nCells + NUM_GROW; i++) {
        for (int j = NUM_GROW; j < y_nCells + NUM_GROW; j++) {
            double x = SimMesh.getx(i);
            double y = SimMesh.gety(j);

//            u[i][j].rho = 1;
//            u[i][j].momentum_x = 1;
//            u[i][j].momentum_y = 0;
//            u[i][j].momentum_z = 0;
//            u[i][j].energy =1;
//            u[i][j].B_x=0;
//            u[i][j].B_y=0;
//            u[i][j].B_z=0;

//            u[i][j].rho = (x+y<2?1.0:0.125);
//            u[i][j].momentum_x = 0;
//            u[i][j].momentum_y = 0;
//            u[i][j].momentum_z = 0;
//            u[i][j].energy = (y+x<2?2.5:0.25);
//            u[i][j].B_x=0;
//            u[i][j].B_y=0;
//            u[i][j].B_z=0;

//            u[i][j].rho = (a_input*x+b_input*y+c_input>0?1.0:0.125);
//            u[i][j].momentum_x = 0;
//            u[i][j].momentum_y = 0;
//            u[i][j].momentum_z = 0;
//            u[i][j].energy = (a_input*x+b_input*y+c_input>0?2.5:0.25);
//            u[i][j].B_x=0;
//            u[i][j].B_y=0;
//            u[i][j].B_z=0;

//            u[i][j].rho = (a_input*x+b_input*y+c_input>0?1.0:0.125);
//            u[i][j].momentum_x = 0;
//            u[i][j].momentum_y = 0;
//            u[i][j].momentum_z = 0;
//            u[i][j].pressure=(a_input*x+b_input*y+c_input>0?1.0:0.1);
//            double shockTube_Bx=0.75,shockTube_By=(a_input*x+b_input*y+c_input>0?1.0:-1.0);
//            u[i][j].B_x=shockTube_Bx*cos(alpha)+shockTube_By* cos(alpha+0.5*M_PI);
//            u[i][j].B_y=shockTube_Bx* sin(alpha)+shockTube_By* sin(alpha+0.5*M_PI);
//            u[i][j].B_z=0;
//            u[i][j].energy = u[i][j].pressure/(C_gamma-1)
//                             +0.5*(u[i][j].momentum_x*u[i][j].momentum_x+u[i][j].momentum_y*u[i][j].momentum_y+u[i][j].momentum_z*u[i][j].momentum_z)/u[i][j].rho
//                             +0.5*(u[i][j].B_x*u[i][j].B_x+u[i][j].B_y*u[i][j].B_y+u[i][j].B_z*u[i][j].B_z);

//            double alpha=0;
//            double a_input=1.0,b_input=0,c_input=-1.0;
//            u[i][j].rho = (a_input*x+b_input*y+c_input<0?1.0:0.125);
//            u[i][j].momentum_x = 0;
//            u[i][j].momentum_y = 0;
//            u[i][j].momentum_z = 0;
//            u[i][j].pressure=(a_input*x+b_input*y+c_input<0?1.0:0.1);
//            double shockTube_Bx=0.75,shockTube_By=(a_input*x+b_input*y+c_input<0?1.0:-1.0);
//            u[i][j].B_x=shockTube_Bx*cos(alpha)+shockTube_By* cos(alpha+0.5*M_PI);
//            u[i][j].B_y=shockTube_Bx* sin(alpha)+shockTube_By* sin(alpha+0.5*M_PI);
//            u[i][j].B_z=0;
//            u[i][j].energy =u[i][j].pressure/(C_gamma-1)
//                             +0.5*(u[i][j].momentum_x*u[i][j].momentum_x+u[i][j].momentum_y*u[i][j].momentum_y+u[i][j].momentum_z*u[i][j].momentum_z)/u[i][j].rho
//                             +0.5*(u[i][j].B_x*u[i][j].B_x+u[i][j].B_y*u[i][j].B_y+u[i][j].B_z*u[i][j].B_z);

            u[i][j].rho = ((x+y)<2?1.0:0.125);
            u[i][j].momentum_x = 0;
            u[i][j].momentum_y = 0;
            u[i][j].momentum_z = 0;
            u[i][j].energy = ((x+y)<2?1.0:0.1)+0.78125;
            u[i][j].B_x=0.53033-0.70710*((x+y)<2?1.0:-1.0);
            u[i][j].B_y=0.53033+0.70710*((x+y)<2?1.0:-1.0);
            u[i][j].B_z=0;
        }
    }
    uPlus1=u;
    rigidBodyBC(u);
//    transmissive(u,NUM_GROW);
    diagonal_transmissive(u,NUM_GROW);
    double t = tStart;
// 初始化结束

    while (t < tStop) {
        dt = ComputeTimeStep(u, dx, dy, 0.4,NUM_GROW);
        if (t + dt > tStop) {
            dt = tStop - t;
            t = tStop;
//        } else if (t<0.5&&(t+dt)>=0.5){
//            dt=0.5-t;
//            t=0.5;
        } else {
            t = t + dt;
        }


        //beginning of y update
        getFlux(flux, u, dy, dt/2,'y',NUM_GROW);
        for (int i = 0; i < x_nCells; i++) {
            for (int j = 0; j < y_nCells; j++) {
                uPlus1[NUM_GROW+i][NUM_GROW+j] = u[NUM_GROW+i][NUM_GROW+j] - (dt/2 / dy) * (flux[i][j+1] - flux[i][j]);
            }
        }
        u=uPlus1;
        rigidBodyBC(u);
//        transmissive(u,NUM_GROW);
        diagonal_transmissive(u,NUM_GROW);
        //end of y update


        //x update
        getFlux(flux, u, dx, dt, 'x',NUM_GROW);
        for (int i = 0; i < x_nCells; i++) {
            for (int j = 0; j < y_nCells; j++) {
                uPlus1[NUM_GROW+i][NUM_GROW+j] = u[NUM_GROW+i][NUM_GROW+j] - (dt / dx) * (flux[i+1][j] - flux[i][j]);
            }
        }
        u=uPlus1;
        rigidBodyBC(u);
//        transmissive(u,NUM_GROW);
        diagonal_transmissive(u,NUM_GROW);
        //end of x update


        //beginning of y update
        getFlux(flux, u, dy, dt/2,'y',NUM_GROW);

        for (int i = 0; i < x_nCells; i++) {
            for (int j = 0; j < y_nCells; j++) {
                uPlus1[NUM_GROW+i][NUM_GROW+j] = u[NUM_GROW+i][NUM_GROW+j] - (dt/2 / dy) * (flux[i][j+1] - flux[i][j]);
            }
        }
        u=uPlus1;
        rigidBodyBC(u);
//        transmissive(u,NUM_GROW);
        diagonal_transmissive(u,NUM_GROW);
        //end of y update

        //coordinate source update//TODO:coordinate
//        for (int i = 0; i < x_nCells; i++) {
//            for (int j = 0; j < y_nCells; j++) {
//                double x = x0 + (i + 0.5-NUM_GROW) * dx;
//                uPlus1[NUM_GROW+i][NUM_GROW+j] = u[NUM_GROW+i][NUM_GROW+j] + dt* u[NUM_GROW+i][NUM_GROW+j].CoordinateSource(1,x);
//            }
//        }
//        transmissive(uPlus1,NUM_GROW);
//        u=uPlus1;

//        DivergencdCleaning(psi,u,psi_flux,dx,dy,dt,NUM_GROW);//TODO:div

        cout<<"time ="<<t<<endl;
    }

    ofstream outfile("\\\\wsl.localhost\\Ubuntu-18.04\\CFD_file\\rigidBody\\Euler.dat");
    if (!outfile)
        return -1;

    for (int i = NUM_GROW; i < x_nCells + NUM_GROW; i++) {
        for (int j = NUM_GROW; j < y_nCells + NUM_GROW; j++) {
            double x = SimMesh.getx(i);
            double y = SimMesh.gety(j);
//            if (get<0>(Phi[i][j])<0) {
//                Equ t(1.0e+5,1.0e+10,1.0e+10,1.0e+10,1.0e+15,1.0e+5,1.0e+5,1.0e+5);
//                outfile << x << " " << y << " " <<t
//                        <<" "<<get<0>(Phi[i][j])
//                        << endl;
//                continue;
//            }
            outfile << x << " " << y << " " << u[i][j]
            <<" "<<get<0>(Phi[i][j])
            << endl;
        }
        outfile<<endl;
    }
    outfile.close();



    ofstream outfile2("\\\\wsl.localhost\\Ubuntu-18.04\\CFD_file\\rigidBody\\Euler_lineout.dat");
    if (!outfile2)
        return -1;

    int resolution=1000;
    double output_dx=(outputEnd_x-outputBegin_x)/(double)resolution;
    double output_dy=(outputEnd_y-outputBegin_y)/(double)resolution;
    for (int i = 0; i < resolution; i++) {
        double x = outputBegin_x+i*output_dx;
        double y = outputBegin_y+i*output_dy;
        double cell_x=(x-x0)/dx-0.5+NUM_GROW;
        double cell_y=(y-y0)/dy-0.5+NUM_GROW;
        int x_1=static_cast<int>(cell_x),y_1=static_cast<int>(cell_y);
        int x_2=x_1+1,y_2=y_1+1;
        Equ result= bilinearInterpolation(u[x_1][y_1],u[x_2][y_1],u[x_1][y_2],u[x_2][y_2],(cell_x-x_1),(cell_y-y_1));
            outfile2 <<(1.5*i/resolution) << " " << result
                    << endl;
    }
    outfile2.close();


//    int tx=10,ty=20;
//    cout<<SimMesh.getx(tx)<<' '<<SimMesh.gety(ty)<<endl;
//    cout<<u[tx][ty][0]<<endl;
//    cout<<get<0>(Phi[tx][ty])<<" the normal is "<<get<1>(Phi[tx][ty]);
//    cout<<"The "<<28<<","<<68<<' '<<"point ";
//    if (  get<0>(Phi[28][68])<=0  &&  get<0>(Phi[28][68])>=(-max(dx,dy)*(NUM_GROW+1))  )
//        for (int k=0;k<8;k++){
//            double T=0.0;
//            cout<<k<<" var :";
//            for (int l = 0; l < Weight[28][68][k].size(); ++l) {
//                tuple<double,double*> t=Weight[28][68][k][l];
//                T=T+get<0>(t)*(*get<1>(t));
//                cout<<"the "<<l<<" is "<<get<0>(t)<<"of"<<*get<1>(t)<<" ";
//            }
//            cout<<endl;
//        }
//    cout<<endl;
    return 0;
}
