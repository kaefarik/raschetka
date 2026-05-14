#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>
#include <windows.h>

using namespace std;

int main() {
    SetConsoleOutputCP(65001); // UTF-8
    SetConsoleCP(65001);

    vector<double> x = {1.415, 1.420, 1.425, 1.430, 1.435, 1.440};
    vector<double> y = {0.888551, 0.889559, 0.890637, 0.891667, 0.892687, 0.893698};
    double xStar = 1.4179;

    double h = 0.005;

    double x2 = 1.420;

    cout << fixed << setprecision(4);

    // ── Производные в точке x2 = 1.376 (индекс 1) ──
    // Центральные разности: y'  = (y[i+1] - y[i-1]) / (2h)
    //                       y'' = (y[i+1] - 2y[i] + y[i-1]) / h^2
    int i = 1;
    double yPrime2  = (y[i + 1] - y[i - 1]) / (2.0 * h);
    double yDPrime2 = (y[i + 1] - 2.0 * y[i] + y[i - 1]) / (h * h);

    cout << "В точке x2 = " << x2 << ":\n";
    cout << "  y'  = " << yPrime2  << "\n";
    cout << "  y'' = " << yDPrime2 << "\n";

    // ── Производные в точке x* = 1.4179 ──
    // x* лежит между x[0]=1.415 и x[1]=1.420
    // Используем формулу Ньютона вперёд от узла x[0]:
    double t = (xStar - x[0]) / h;
    double dy1 = y[1] - y[0];
    double dy2 = y[2] - y[1];
    double d2y1 = dy2 - dy1;

    double yPrimeStar = (dy1 + (2.0 * t - 1.0) / 2.0 * d2y1) / h;
    double yDPrimeStar = d2y1 / (h * h);

    cout << "\nВ точке x* = " << xStar << ":\n";
    cout << "  y'  = " << yPrimeStar  << "\n";
    cout << "  y'' = " << yDPrimeStar << "\n";

    return 0;
}
