#include <iostream>
#include <vector>
#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace boost::numeric::odeint;

// 定义自定义类型 Equ
class Equ {
public:
    double value;

    Equ() : value(0.0) {}
    Equ(double v) : value(v) {}

    // 支持加法
    Equ operator+(const Equ &rhs) const {
        return Equ(this->value + rhs.value);
    }

    // 支持减法
    Equ operator-(const Equ &rhs) const {
        return Equ(this->value - rhs.value);
    }

    // 支持乘法
    Equ operator*(double rhs) const {
        return Equ(this->value * rhs);
    }

    // 支持除法
    Equ operator/(double rhs) const {
        return Equ(this->value / rhs);
    }

    // 支持赋值
    Equ& operator=(const Equ &rhs) {
        if (this != &rhs) {
            this->value = rhs.value;
        }
        return *this;
    }

    // 支持累加
    Equ& operator+=(const Equ &rhs) {
        this->value += rhs.value;
        return *this;
    }

    // 支持累减
    Equ& operator-=(const Equ &rhs) {
        this->value -= rhs.value;
        return *this;
    }

    // 支持乘以标量
    Equ& operator*=(double rhs) {
        this->value *= rhs;
        return *this;
    }

    // 支持除以标量
    Equ& operator/=(double rhs) {
        this->value /= rhs;
        return *this;
    }

    friend ostream& operator<<(ostream &os, const Equ &equ) {
        os << equ.value;
        return os;
    }
};

// 定义状态类型
typedef vector<Equ> state_type;

// 定义微分方程
void dydt(const state_type &y, state_type &dy, double t) {
    dy[0] = -2 * y[0];
}

int main() {
    // 初始条件
    state_type y(1);
    y[0] = Equ(1.0);

    // 时间范围和步长
    double t = 0.0;
    double t_end = 5.0;
    double dt = 0.1;

    // 结果存储
    vector<double> times;
    vector<state_type> values;

    // 创建 RK4 步进器
    runge_kutta4<state_type> stepper;

    // 时间积分
    while (t < t_end) {
        times.push_back(t);
        values.push_back(y);
        stepper.do_step(dydt, y, t, dt);
        t += dt;
    }

    // 输出结果
    for (size_t i = 0; i < times.size(); ++i) {
        cout << "t: " << times[i] << ", y: " << values[i][0] << endl;
    }

    return 0;
}
