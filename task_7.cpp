#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <string>
#include <windows.h>

using namespace std;

double sumPower(const vector<double>& arr, int power) {
    double s = 0;
    for (double x : arr) s += pow(x, power);
    return s;
}

double sumXYPower(const vector<double>& x, const vector<double>& y, int xPower) {
    double s = 0;
    for (int i = 0; i < (int)x.size(); i++)
        s += pow(x[i], xPower) * y[i];
    return s;
}

// Метод Гаусса для системы n×n
vector<double> gaussElim(vector<vector<double>> A, vector<double> B) {
    int n = B.size();
    for (int i = 0; i < n; i++) {
        double pivot = A[i][i];
        for (int j = i; j < n; j++) A[i][j] /= pivot;
        B[i] /= pivot;
        for (int k = i + 1; k < n; k++) {
            double f = A[k][i];
            for (int j = i; j < n; j++) A[k][j] -= f * A[i][j];
            B[k] -= f * B[i];
        }
    }
    vector<double> X(n, 0);
    for (int i = n - 1; i >= 0; i--) {
        X[i] = B[i];
        for (int j = i + 1; j < n; j++) X[i] -= A[i][j] * X[j];
    }
    return X;
}

double rSquared(const vector<double>& yTrue, const vector<double>& yPred) {
    double yMean = 0;
    for (double v : yTrue) yMean += v;
    yMean /= yTrue.size();
    double ssRes = 0, ssTot = 0;
    for (int i = 0; i < (int)yTrue.size(); i++) {
        ssRes += (yTrue[i] - yPred[i]) * (yTrue[i] - yPred[i]);
        ssTot += (yTrue[i] - yMean)    * (yTrue[i] - yMean);
    }
    return 1.0 - ssRes / ssTot;
}

int main() {
    SetConsoleOutputCP(65001);   // UTF-8
    SetConsoleCP(65001);
    vector<double> x = {0, 0.12, 0.19, 0.35, 0.4, 0.45, 0.62, 0.71, 0.84, 0.91, 1.0};
    vector<double> y = {-1.0, -1.0, -0.9, -0.5, -0.7, -0.6, -0.3, -0.5, 0.4, 0.8, 1.2};
    int n = x.size();

    cout << fixed;

    // ── Линейная аппроксимация ──
    double Sx  = sumPower(x, 1);
    double Sy  = sumXYPower(x, y, 0);  // = sum(y)
    double Sx2 = sumPower(x, 2);
    double Sxy = sumXYPower(x, y, 1);

    // Пересчитаем Sy как просто сумму y
    Sy = 0; for (double v : y) Sy += v;

    cout << "=== ЛИНЕЙНАЯ АППРОКСИМАЦИЯ ===\n";
    cout << setprecision(4);
    cout << "Σx  = " << Sx  << "\n";
    cout << "Σy  = " << Sy  << "\n";
    cout << "Σx2 = " << Sx2 << "\n";
    cout << "Σxy = " << Sxy << "\n";

    double det   = n * Sx2 - Sx * Sx;
    double aLin  = (Sy * Sx2 - Sx * Sxy) / det;
    double bLin  = (n * Sxy - Sx * Sy)   / det;

    cout << "\nСистема уравнений:\n";
    cout << n << "a + " << Sx << "b = " << Sy << "\n";
    cout << Sx << "a + " << Sx2 << "b = " << Sxy << "\n";
    cout << "\nРешение: a = " << aLin << ", b = " << bLin << "\n";
    cout << "Уравнение: y = " << aLin << " + " << bLin << "x\n";

    // ── Квадратичная аппроксимация ──
    double Sx3  = sumPower(x, 3);
    double Sx4  = sumPower(x, 4);
    double Sx2y = sumXYPower(x, y, 2);

    cout << "\n=== КВАДРАТИЧНАЯ АППРОКСИМАЦИЯ ===\n";
    cout << setprecision(6);
    cout << "Σx3  = " << Sx3  << "\n";
    cout << "Σx4  = " << Sx4  << "\n";
    cout << "Σx2y = " << Sx2y << "\n";

    cout << setprecision(4);
    cout << "\nСистема уравнений:\n";
    cout << n << "a + " << Sx << "b + " << Sx2 << "c = " << Sy << "\n";
    cout << Sx << "a + " << Sx2 << "b + " << Sx3 << "c = " << Sxy << "\n";
    cout << Sx2 << "a + " << Sx3 << "b + " << Sx4 << "c = " << Sx2y << "\n";

    vector<vector<double>> A = {
        {(double)n, Sx,  Sx2},
        {Sx,  Sx2, Sx3},
        {Sx2, Sx3, Sx4}
    };
    vector<double> B = {Sy, Sxy, Sx2y};
    auto sol = gaussElim(A, B);
    double aQ = sol[0], bQ = sol[1], cQ = sol[2];

    cout << "\nРешение: a = " << aQ << ", b = " << bQ << ", c = " << cQ << "\n";
    cout << "Уравнение: y = " << aQ << " + " << bQ << "x + " << cQ << "x2\n";

    // ── Оценка качества ──
    vector<double> yLin(n), yQuad(n);
    for (int i = 0; i < n; i++) {
        yLin[i]  = aLin + bLin * x[i];
        yQuad[i] = aQ   + bQ   * x[i] + cQ * x[i] * x[i];
    }

    cout << setprecision(6);
    cout << "\n=== ОЦЕНКА КАЧЕСТВА ===\n";
    cout << "R2 для линейной модели:     " << rSquared(y, yLin)  << "\n";
    cout << "R2 для квадратичной модели: " << rSquared(y, yQuad) << "\n";

    // ── Сравнение значений ──
    cout << "\n=== СРАВНЕНИЕ ЗНАЧЕНИЙ ===\n";
    cout << setw(8)  << "x"
         << setw(8)  << "y"
         << setw(10) << "y_лин"
         << setw(10) << "y_квад" << "\n";
    cout << string(40, '-') << "\n";
    for (int i = 0; i < n; i++) {
        cout << setprecision(2) << setw(8) << x[i]
             << setw(8) << y[i]
             << setprecision(4)
             << setw(10) << yLin[i]
             << setw(10) << yQuad[i] << "\n";
    }

    return 0;
}
