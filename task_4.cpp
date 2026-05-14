#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <string>
#include <limits>
#include <windows.h>

using namespace std;

const double EPS = 0.01;
const double INF = numeric_limits<double>::infinity();

// f(x) = x + lg(x) - 0.5
double f(double x) {
    if (x <= 0) return NAN;
    return x + log10(x) - 0.5;
}

// g(x) = 0.5 - lg(x)  (итерационная функция: x = g(x))
double g(double x) {
    if (x <= 0) return NAN;
    return 0.5 - log10(x);
}

double derivative(double (*func)(double), double x, double h = 1e-7) {
    return (func(x + h) - func(x - h)) / (2.0 * h);
}

// ─── Метод половинного деления ───────────────────────────────────────────────

struct HalfRow { int n; double a, b, xn, fxn; };

pair<double, vector<HalfRow>> halfDivision(double e) {
    double a = 0.5;
    double b_ = 1;

    vector<HalfRow> rows;
    int n = (int)ceil(log2((b_ - a) / (2.0 * e)));
    for (int i = 0; i < n; i++) {
        double xn = (a + b_) / 2.0;
        double fxn = f(xn);
        rows.push_back({i + 1,
            round(a   * 1e6) / 1e6,
            round(b_  * 1e6) / 1e6,
            round(xn  * 1e6) / 1e6,
            round(fxn * 1e6) / 1e6});
        if (f(a) * fxn < 0) b_ = xn;
        else                 a  = xn;
    }
    return {(a + b_) / 2.0, rows};
}

// ─── Метод простой итерации ──────────────────────────────────────────────────

struct IterRow2 { int n; double x0, diff; };

pair<double, vector<IterRow2>> simpleIteration(double e) {
    double x0 = -1, q = -1;
    double step = 0.1;
    for (double i = 0.1; i < 10.0 && x0 < 0; i = round((i + step) * 100) / 100) {
        for (double k = round((i + step) * 100) / 100; k < 10.0 && x0 < 0;
             k = round((k + step) * 100) / 100) {
            double fi = f(i), fk = f(k);
            if (!isnan(fi) && !isnan(fk) && fi * fk < 0) {
                double dgi = fabs(derivative(g, i));
                double dgk = fabs(derivative(g, k));
                if (dgi < 1.0 && dgk < 1.0) {
                    q  = max(dgi, dgk);
                    x0 = (i + k) / 2.0;
                }
            }
        }
    }

    vector<IterRow2> rows;
    double x1 = g(x0);
    int n = 1;
    while (fabs(x1 - x0) > ((1.0 - q) / q) * e) {
        rows.push_back({n, round(x0 * 1e6) / 1e6, round(fabs(x1 - x0) * 1e6) / 1e6});
        x0 = x1;
        x1 = g(x0);
        n++;
    }
    rows.push_back({n, round(x0 * 1e6) / 1e6, round(fabs(x1 - x0) * 1e6) / 1e6});
    return {x1, rows};
}

// ─── Метод Ньютона ───────────────────────────────────────────────────────────

struct NewtonRow { int n; double x0, fn, dfn, x1; };

pair<double, vector<NewtonRow>> newtonMethod(double e) {
    double x0 = -1;
    double step = 0.1;
    for (double i = 0.1; i < 10.0 && x0 < 0; i = round((i + step) * 100) / 100) {
        for (double k = round((i + step) * 100) / 100; k < 10.0 && x0 < 0;
             k = round((k + step) * 100) / 100) {
            double fi = f(i), fk = f(k);
            if (!isnan(fi) && !isnan(fk) && fi * fk < 0)
                x0 = (i + k) / 2.0;
        }
    }

    vector<NewtonRow> rows;
    double x1 = x0 - f(x0) / derivative(f, x0);
    int n = 1;
    while (fabs(x1 - x0) > e) {
        double fn = f(x0), dfn = derivative(f, x0);
        rows.push_back({n,
            round(x0  * 1e6) / 1e6,
            round(fn  * 1e6) / 1e6,
            round(dfn * 1e6) / 1e6,
            round(x1  * 1e6) / 1e6});
        x0 = x1;
        x1 = x0 - f(x0) / derivative(f, x0);
        n++;
    }
    rows.push_back({n,
        round(x0 * 1e6) / 1e6,
        round(f(x0) * 1e6) / 1e6,
        round(derivative(f, x0) * 1e6) / 1e6,
        round(x1 * 1e6) / 1e6});
    return {x1, rows};
}

// ─── main ────────────────────────────────────────────────────────────────────

int main() {
    SetConsoleOutputCP(65001);   // UTF-8
    SetConsoleCP(65001);
    cout << fixed << setprecision(6);

    auto [rootHalf, rowsHalf] = halfDivision(EPS);
    auto [rootIter, rowsIter] = simpleIteration(EPS);
    auto [rootNewt, rowsNewt] = newtonMethod(EPS);

    cout << "Простая итерация:   " << rootIter << "\n";
    cout << "Половинное деление: " << rootHalf << "\n";
    cout << "Метод Ньютона:      " << rootNewt << "\n\n";

    // Метод половинного деления
    cout << string(60, '-') << "\nМетод половинного деления\n" << string(60, '-') << "\n";
    cout << setw(5) << "n"
         << setw(12) << "a"
         << setw(12) << "b"
         << setw(12) << "xn"
         << setw(12) << "f(xn)" << "\n";
    for (auto& r : rowsHalf)
        cout << setw(5) << r.n
             << setw(12) << r.a
             << setw(12) << r.b
             << setw(12) << r.xn
             << setw(12) << r.fxn << "\n";
    cout << "Корень: x* ≈ " << rootHalf << "\n\n";

    // Метод простой итерации
    cout << string(60, '-') << "\nМетод простой итерации\n" << string(60, '-') << "\n";
    cout << setw(5) << "n"
         << setw(12) << "xn"
         << setw(16) << "|xn - xn-1|" << "\n";
    for (auto& r : rowsIter)
        cout << setw(5) << r.n
             << setw(12) << r.x0
             << setw(16) << r.diff << "\n";
    cout << "Корень: x* ≈ " << rootIter << "\n\n";

    // Метод Ньютона
    cout << string(60, '-') << "\nМетод Ньютона\n" << string(60, '-') << "\n";
    cout << setw(5) << "n"
         << setw(12) << "xn-1"
         << setw(12) << "f(xn-1)"
         << setw(12) << "f'(xn-1)"
         << setw(12) << "xn" << "\n";
    for (auto& r : rowsNewt)
        cout << setw(5) << r.n
             << setw(12) << r.x0
             << setw(12) << r.fn
             << setw(12) << r.dfn
             << setw(12) << r.x1 << "\n";
    cout << "Корень: x* ≈ " << rootNewt << "\n";

    return 0;
}
