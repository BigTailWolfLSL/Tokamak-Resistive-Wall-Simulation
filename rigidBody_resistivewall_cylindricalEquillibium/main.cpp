#include "bili_interpolation.h"
#include <iostream>

using namespace std;

int main() {
    double a=1,b=2,c=3,d=4;
    double result= bilinearInterpolation(a,b,c,d,0.5,0);
    cout<<result<<endl;
    return 0;
}
