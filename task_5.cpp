#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <string>
#include <windows.h>

using namespace std;

// ─── Интерполяция Лагранжа ───────────────────────────────────────────────────

double lagrange(double x, const vector<double>& xs, const vector<double>& ys) {
    int n = xs.size();
    double result = 0.0;
    for (int i = 0; i < n; i++) {
        double Li = ys[i];
        for (int j = 0; j < n; j++) {
            if (j != i)
                Li *= (x - xs[j]) / (xs[i] - xs[j]);
        }
        result += Li;
    }
    return result;
}

// ─── Разделённые разности Ньютона ────────────────────────────────────────────

vector<double> newtonCoefficients(const vector<double>& xs, vector<double> ys) {
    int n = xs.size();
    for (int j = 1; j < n; j++)
        for (int i = n - 1; i >= j; i--)
            ys[i] = (ys[i] - ys[i - 1]) / (xs[i] - xs[i - j]);
    return ys;
}

double newtonInterpolation(double x, const vector<double>& xs, const vector<double>& coef) {
    int n = xs.size();
    double result = coef[0];
    double product = 1.0;
    for (int i = 1; i < n; i++) {
        product *= (x - xs[i - 1]);
        result  += coef[i] * product;
    }
    return result;
}

double errorEstimate(double x, const vector<double>& xs, const vector<double>& coef) {
    int n = xs.size();
    double product = 1.0;
    for (double xi : xs) product *= (x - xi);
    return coef[n - 1] * product;
}

// ─── main ────────────────────────────────────────────────────────────────────

int main() {
    SetConsoleOutputCP(65001);   // UTF-8
    SetConsoleCP(65001);
    vector<double> xData = {1.415, 1.420, 1.425, 1.430, 1.435, 1.440};
    vector<double> yData = {0.888551, 0.889559, 0.890637, 0.891667, 0.892687, 0.893698};
    double xStar = 1.4179;

    auto newtonCoef = newtonCoefficients(xData, yData);

    double yLagrange = lagrange(xStar, xData, yData);
    double yNewton   = newtonInterpolation(xStar, xData, newtonCoef);
    double errNewton = errorEstimate(xStar, xData, newtonCoef);

    cout << fixed;
    cout << string(60, '=') << "\n";
    cout << "ИНТЕРПОЛЯЦИЯ ПО УЗЛАМ\n";
    cout << string(60, '=') << "\n";

    cout << "Узлы x_i:    ";
    for (double x : xData) cout << setprecision(3) << x << "  ";
    cout << "\n";

    cout << "Значения y_i: ";
    for (double y : yData) cout << setprecision(5) << y << "  ";
    cout << "\n";

    cout << "Точка x* = " << setprecision(4) << xStar << "\n";
    cout << string(60, '-') << "\n";

    cout << setprecision(8);
    cout << "Значение по Лагранжу:   y ≈ " << yLagrange << "\n";
    cout << "Значение по Ньютону:    y ≈ " << yNewton   << "\n";
    cout << string(60, '-') << "\n";
    cout << "Оценка погрешности (Ньютон): Δy ≈ " << scientific << setprecision(2) << errNewton << "\n";
    cout << string(60, '=') << "\n\n";

    // Коэффициенты многочлена Ньютона
    cout << "Коэффициенты многочлена Ньютона (разделённые разности):\n";
    string labels[] = {"f[x0]        ", "f[x0,x1]     ", "f[x0,x1,x2]  ",
                        "f[x0,...,x3] ", "f[x0,...,x4] ", "f[x0,...,x5] "};
    cout << fixed << setprecision(8);
    for (int i = 0; i < (int)newtonCoef.size(); i++)
        cout << labels[i] << "= " << newtonCoef[i] << "\n";

    cout << "\nМногочлен Ньютона:\n";
    cout << "N(x) = " << newtonCoef[0] << "\n";
    for (int i = 1; i < (int)newtonCoef.size(); i++) {
        cout << "     + " << newtonCoef[i];
        for (int j = 0; j < i; j++)
            cout << "·(x - " << setprecision(3) << xData[j] << ")";
        cout << setprecision(8) << "\n";
    }

    return 0;
}
