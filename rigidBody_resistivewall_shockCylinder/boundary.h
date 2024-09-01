#ifndef boundary_H
#define boundary_H

/**
 * Boundary Condition function in dim s.
 */

#include <vector>
using namespace std;


/**
 * 1D boundary condition
 */
template<class T>
void periodic(vector<T> &uPlus1) {
    //periodic boundary condition
    uPlus1[0] = uPlus1[uPlus1.size() - 2];
    uPlus1[uPlus1.size() - 1] = uPlus1[1];
}

template<class T>
void transmissive(vector<T> &uPlus1) {
    //transmissive boundary condition
    uPlus1[0] = uPlus1[1];
    uPlus1[uPlus1.size() - 1] = uPlus1[uPlus1.size() - 2];
}

template<class T>
void reflective(vector<T> &uPlus1) {
    //reflective boundary condition
    uPlus1[0] = -uPlus1[1];
    uPlus1[uPlus1.size() - 1] = -uPlus1[uPlus1.size() - 2];
}

template<class T>
void fixed(vector<T> &uPlus1, double value = 0) {
    //reflective boundary condition
    uPlus1[0] = value;
    uPlus1[uPlus1.size() - 1] = value;
}

/**
 * 2D boundary condition
 */

template<class T>
void periodic(vector<vector<T>> &uPlus1, int NUM_GROW) {
    int nx = uPlus1.size();
    int ny = uPlus1[0].size();

    // 处理x方向的周期性边界
    for (int l = 0; l < NUM_GROW; l++) {
        for (int j = 0; j < ny; j++) {
            uPlus1[nx - NUM_GROW + l][j] = uPlus1[NUM_GROW + l][j];// 从上到下
            uPlus1[l][j] = uPlus1[nx - 2 * NUM_GROW + l][j];       // 从下到上
        }
    }

    // 处理y方向的周期性边界
    for (int l = 0; l < NUM_GROW; l++) {
        for (int i = 0; i < nx; i++) {
            uPlus1[i][ny - NUM_GROW + l] = uPlus1[i][NUM_GROW + l];// 从左到右
            uPlus1[i][l] = uPlus1[i][ny - 2 * NUM_GROW + l];       // 从右到左
        }
    }
}

template<class T>
void transmissive(vector<vector<T>> &uPlus1, int NUM_GROW) {
    int nx = uPlus1.size();
    int ny = uPlus1[0].size();// 应该使用uPlus1[0].size()获取第一行的长度，以确保是y维度的大小

    // 处理y方向的边界
    for (int l = 0; l < NUM_GROW; l++) {
        for (int i = NUM_GROW; i < nx - NUM_GROW; i++) {
            uPlus1[i][NUM_GROW - l - 1] = uPlus1[i][NUM_GROW - l];          // 左边界
            uPlus1[i][ny - NUM_GROW + l] = uPlus1[i][ny - NUM_GROW + l - 1];// 右边界
        }
    }

    // 处理x方向的边界
    for (int l = 0; l < NUM_GROW; l++) {
        for (int j = 0; j < ny; j++) {
            uPlus1[NUM_GROW - l - 1][j] = uPlus1[NUM_GROW - l][j];          // 上边界
            uPlus1[nx - NUM_GROW + l][j] = uPlus1[nx - NUM_GROW + l - 1][j];// 下边界
        }
    }
}

/**
 * Kelvin-Helmholtz Instability的边界条件，reflective on top and bottom, periodic on left and right
 * @tparam T 肯定是个Equ类型，这里特定是MHD的8个量
 * @param uPlus1
 * @param NUM_GROW
 */
template<class T>
void customBoundary(vector<vector<T>> &uPlus1, int NUM_GROW) {
    int nx = uPlus1.size();
    int ny = uPlus1[0].size();

    // 处理x方向的周期性边界条件
    for (int j = 0; j < ny; j++) {
        for (int l = 0; l < NUM_GROW; l++) {
            // 左边界复制到右边界
            uPlus1[nx - NUM_GROW + l][j] = uPlus1[NUM_GROW + l][j];
            // 右边界复制到左边界
            uPlus1[l][j] = uPlus1[nx - 2 * NUM_GROW + l][j];
        }
    }

    // 处理y方向的反射性边界条件
    for (int i = 0; i < nx; i++) {
        // 反射性边界条件，y上端
        for (int l = 0; l < NUM_GROW; l++) {
            uPlus1[i][NUM_GROW - l - 1] = uPlus1[i][NUM_GROW + l];
            // 反射边界对速度vy取反
            uPlus1[i][NUM_GROW - l - 1][2] *= -1;// 反转vy
            uPlus1[i][NUM_GROW - l - 1].B_y *= -1;

            uPlus1[i][ny - NUM_GROW + l] = uPlus1[i][ny - NUM_GROW - 1 - l];
            // 反射边界对速度vy取反
            uPlus1[i][ny - NUM_GROW + l][2] *= -1;// 反转vy
            uPlus1[i][ny - NUM_GROW + l].B_y *= -1;
            //磁场不能取反，否则创造磁单极
        }
    }
}
template<class T>
void customBoundarypsi(vector<vector<T>> &uPlus1, int NUM_GROW) {
    int nx = uPlus1.size();
    int ny = uPlus1[0].size();

    // 处理x方向的周期性边界条件
    for (int j = 0; j < ny; j++) {
        for (int l = 0; l < NUM_GROW; l++) {
            // 左边界复制到右边界
            uPlus1[nx - NUM_GROW + l][j] = uPlus1[NUM_GROW + l][j];
            // 右边界复制到左边界
            uPlus1[l][j] = uPlus1[nx - 2 * NUM_GROW + l][j];
        }
    }

    // 处理y方向的反射性边界条件
    for (int i = 0; i < nx; i++) {
        // 反射性边界条件，y上端
        for (int l = 0; l < NUM_GROW; l++) {
            uPlus1[i][NUM_GROW - l - 1] = uPlus1[i][NUM_GROW + l];
            // 反射边界对速度vy取反
            //            uPlus1[i][NUM_GROW - l - 1][2] *= -1; // 反转vy
            //            uPlus1[i][NUM_GROW - l - 1].B_y *= -1;

            uPlus1[i][ny - NUM_GROW + l] = uPlus1[i][ny - NUM_GROW - 1 - l];
            // 反射边界对速度vy取反
            //            uPlus1[i][ny - NUM_GROW + l][2] *= -1; // 反转vy
            //            uPlus1[i][ny - NUM_GROW + l].B_y *= -1;
            //磁场不能取反，否则创造磁单极
        }
    }
}
template<class T>
void customBoundary2(vector<vector<T>> &uPlus1, int NUM_GROW) {
    int nx = uPlus1.size();
    int ny = uPlus1[0].size();

    // 处理x方向的反射性边界条件
    for (int j = 0; j < ny; j++) {
        for (int l = 0; l < NUM_GROW; l++) {
            // 左边界反射到右边界
            uPlus1[nx - NUM_GROW + l][j] = uPlus1[nx - NUM_GROW - 1 - l][j];
            // 右边界反射到左边界
            uPlus1[l][j] = uPlus1[2 * NUM_GROW - 1 - l][j];
        }
    }

    // 处理y方向的周期性边界条件
    for (int i = 0; i < nx; i++) {
        for (int l = 0; l < NUM_GROW; l++) {
            // 上边界复制到下边界
            uPlus1[i][ny - NUM_GROW + l] = uPlus1[i][NUM_GROW + l];
            // 下边界复制到上边界
            uPlus1[i][l] = uPlus1[i][ny - 2 * NUM_GROW + l];
        }
    }
}

/**
 * diagonal下 波速y=x方向的 transmissive条件
 * @tparam T
 * @param uPlus1
 * @param NUM_GROW
 */
template<class T>
void diagonal_transmissive(vector<vector<T>> &uPlus1, int NUM_GROW) {
    int nx = uPlus1.size();
    int ny = uPlus1[0].size();

    for (int l = 0; l < NUM_GROW; l++) {
        for (int i = 0; i < (nx - 2 * NUM_GROW + 2 * l); i++) {
            uPlus1[NUM_GROW - 1 - l + i][ny - NUM_GROW + l] = uPlus1[NUM_GROW - l + i][ny - NUM_GROW + l - 1];
            uPlus1[nx - NUM_GROW + l - i][NUM_GROW - 1 - l] = uPlus1[nx - NUM_GROW + l - 1 - i][NUM_GROW - l];
        }

        for (int i = 0; i < (ny - 2 * NUM_GROW + 2 * l); i++) {
            uPlus1[NUM_GROW - 1 - l][ny - NUM_GROW + l - i] = uPlus1[NUM_GROW - l][ny - NUM_GROW + l - 1 - i];
            uPlus1[nx - NUM_GROW + l][NUM_GROW - 1 - l + i] = uPlus1[nx - NUM_GROW + l - 1][NUM_GROW - l + i];
        }

        uPlus1[NUM_GROW - 1 - l][NUM_GROW - 1 - l] = uPlus1[NUM_GROW - 1 - l + 1][NUM_GROW - 1 - l + 1];
        uPlus1[NUM_GROW - 1 - l + 1][NUM_GROW - 1 - l] = uPlus1[NUM_GROW - 1 - l + 1][NUM_GROW - 1 - l + 1];
        uPlus1[NUM_GROW - 1 - l][NUM_GROW - 1 - l + 1] = uPlus1[NUM_GROW - 1 - l + 1][NUM_GROW - 1 - l + 1];

        uPlus1[nx - NUM_GROW + l][ny - NUM_GROW + l] = uPlus1[nx - NUM_GROW + l - 1][ny - NUM_GROW + l - 1];
        uPlus1[nx - NUM_GROW + l - 1][ny - NUM_GROW + l] = uPlus1[nx - NUM_GROW + l - 1][ny - NUM_GROW + l - 1];
        uPlus1[nx - NUM_GROW + l][ny - NUM_GROW + l - 1] = uPlus1[nx - NUM_GROW + l - 1][ny - NUM_GROW + l - 1];
    }
}

/**
 * MHD reflective boundary for Brio-Wu test on x
 * @tparam T
 * @param uPlus1
 * @param NUM_GROW
 */
template<class T>
void Brio_Wu_reflectiveBC(vector<vector<T>> &uPlus1, int NUM_GROW) {
    int nx = uPlus1.size();
    int ny = uPlus1[0].size();


//    for (int l = 0; l < NUM_GROW; l++) {
//        for (int j = NUM_GROW; j < ny-NUM_GROW; j++) {
//            // 左边界反射到左边界有效部分
//            uPlus1[nx - NUM_GROW + l][j] = uPlus1[nx - NUM_GROW - 1 - l][j];
//            uPlus1[nx - NUM_GROW + l][j].B_x = 1.5-uPlus1[nx - NUM_GROW + l][j].B_x;
//            uPlus1[nx - NUM_GROW + l][j].momentum_x *= -1;
//            // 右边界反射到右边界有效部分
//            uPlus1[NUM_GROW-1-l][j] = uPlus1[NUM_GROW+l][j];
//            uPlus1[NUM_GROW-1-l][j].B_x = 1.5-uPlus1[NUM_GROW-1-l][j].B_x;
//            uPlus1[NUM_GROW-1-l][j].momentum_x *= -1;
//        }
//        for (int j=0;j<nx;j++){//y轴两端更新
//            double C=(j<(nx/2))?-1.0:1.0;
//            uPlus1[j][NUM_GROW-1-l]=uPlus1[j][NUM_GROW+l];
//            uPlus1[j][NUM_GROW-1-l].momentum_y*=-1;
////            uPlus1[j][NUM_GROW-1-l].B_y=2*C-uPlus1[j][NUM_GROW-1-l].B_y;
//            uPlus1[j][ny-NUM_GROW+l]=uPlus1[j][ny-NUM_GROW-1-l];
//            uPlus1[j][ny-NUM_GROW+l].momentum_y*=-1;
////            uPlus1[j][ny-NUM_GROW+l].B_y=2*C-uPlus1[j][ny-NUM_GROW+l].B_y;
//        }
//    }

    for (int l = 0; l < NUM_GROW; l++) {
        for (int i = NUM_GROW; i < nx - NUM_GROW; i++) {
            double C=0;//(i<(nx/2))?1.0:-1.0;
            uPlus1[i][NUM_GROW - l - 1] = uPlus1[i][NUM_GROW + l];          // y下边界
            uPlus1[i][NUM_GROW - l - 1].momentum_y*=-1;
            uPlus1[i][NUM_GROW - l - 1].B_y=2*C-uPlus1[i][NUM_GROW - l - 1].B_y;
            uPlus1[i][ny - NUM_GROW + l] = uPlus1[i][ny - NUM_GROW - l - 1];// y上边界
            uPlus1[i][ny - NUM_GROW + l].momentum_y*=-1;
            uPlus1[i][ny - NUM_GROW + l].B_y=2*C-uPlus1[i][ny - NUM_GROW + l].B_y;
        }
    }
    for (int l = 0; l < NUM_GROW; l++) {
        for (int j = 0; j < ny; j++) {
            uPlus1[NUM_GROW - l - 1][j] = uPlus1[NUM_GROW + l][j];          // x左边界
            uPlus1[NUM_GROW - l - 1][j].momentum_x*=-1;
            uPlus1[NUM_GROW - l - 1][j].B_x=1.5-uPlus1[NUM_GROW - l - 1][j].B_x;
            uPlus1[nx - NUM_GROW + l][j] = uPlus1[nx - NUM_GROW - l - 1][j];// x右边界
            uPlus1[nx - NUM_GROW + l][j].momentum_x*=-1;
            uPlus1[nx - NUM_GROW + l][j].B_x=1.5-uPlus1[nx - NUM_GROW + l][j].B_x;
        }
    }
}

#endif//boundary_H