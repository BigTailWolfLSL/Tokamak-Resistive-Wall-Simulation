#ifndef RIGIDBODY_BILI_INTERPOLATION_H
#define RIGIDBODY_BILI_INTERPOLATION_H

using namespace std;

template<class T>
T linearInterpolation(T& v0,T& v1,double r){
    return r*v1+(1-r)*v0;
}

template<class T>
T bilinearInterpolation(T& vDL,T& vDR, T& vUL, T& vUR,double rx,double ry){
    //x horizontal; y vertical; D/U down or up; L/R left or right
    T vD= linearInterpolation(vDL,vDR,rx);
    T vU= linearInterpolation(vUL,vUR,rx);
    return linearInterpolation(vD,vU,ry);
}

#endif//RIGIDBODY_BILI_INTERPOLATION_H
