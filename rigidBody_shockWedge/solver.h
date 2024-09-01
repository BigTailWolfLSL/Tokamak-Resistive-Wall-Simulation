#ifndef RIGIDBODY_SOLVER_H
#define RIGIDBODY_SOLVER_H

#include <cmath>
#include "MHDvec.h"

using namespace std;

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
//            r=delta_l[4]/(delta_r[4]+1.0e-20);
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
//            r=delta_l[4]/(delta_r[4]+1.0e-20);
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
            r=delta_l[4]/(delta_r[4]+1.0e-20);
            for (int l = 0; l < delta_l.n; l++) {
                //                r = delta_l[l] / (delta_r[l] + 1.0e-20);//(delta_r[l] != 0 ? delta_l[l] / delta_r[l] : 1.0e+20);
                deltaHead_l[l]=0.5*(1+omega)*delta_l[l]+0.5*(1-omega)*delta_r[l];
                deltaHead_l[l]=deltaHead_l[l]* limiter(r);
            }

        }
        {
            Equ delta_l = u[i][j] - u[i][j-1];
            Equ delta_r = u[i][j+1] - u[i][j];
            double r;
            r=delta_l[4]/(delta_r[4]+1.0e-20);
            for (int l = 0; l < delta_l.n; l++) {
                //                r = delta_l[l] / (delta_r[l] + 1.0e-20);//(delta_r[l] != 0 ? delta_l[l] / delta_r[l] : 1.0e+20);
                deltaHead_r[l]=0.5*(1+omega)*delta_l[l]+0.5*(1-omega)*delta_r[l];
                deltaHead_r[l]=deltaHead_r[l]* limiter(r);//VanLeer
            }
        }

        Equ ul=(u[i][j-1]+0.5*deltaHead_l)-0.5*dt/dx*((u[i][j-1]+0.5*deltaHead_l).y_fluxFunc() -(u[i][j-1]-0.5*deltaHead_l).y_fluxFunc());
        Equ ur=(u[i][j]-0.5*deltaHead_r)-0.5*dt/dx*((u[i][j]+0.5*deltaHead_r).y_fluxFunc() -(u[i][j]-0.5*deltaHead_r).y_fluxFunc());
        return HLLC_flux2(ul,ur,dx,dt,direction);
    }
}

#endif//RIGIDBODY_SOLVER_H
