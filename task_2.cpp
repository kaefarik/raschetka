#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <string>
#include <windows.h>

using namespace std;

typedef vector<double> Row;
typedef vector<Row> Matrix;

Row division(Row line, double n) {
    for (auto& v : line) v /= n;
    return line;
}

// ─── Метод прогонки ──────────────────────────────────────────────────────────

vector<double> progonka(const Matrix& slay) {
    int n = slay.size();

    vector<double> a(n, 0), b(n), c(n, 0), d(n);

    for (int i = 0; i < n; i++) {
        b[i] = slay[i][i];
        d[i] = slay[i][n];
        if (i > 0)     a[i] = slay[i][i - 1];
        if (i < n - 1) c[i] = slay[i][i + 1];
    }

    vector<double> alpha(n, 0), beta(n, 0);
    alpha[0] = -c[0] / b[0];
    beta[0]  =  d[0] / b[0];

    for (int i = 1; i < n; i++) {
        double denom = b[i] + a[i] * alpha[i - 1];
        alpha[i] = (i < n - 1) ? -c[i] / denom : 0.0;
        beta[i]  = (d[i] - a[i] * beta[i - 1]) / denom;
    }

    vector<double> x(n);
    x[n - 1] = beta[n - 1];
    for (int i = n - 2; i >= 0; i--)
        x[i] = alpha[i] * x[i + 1] + beta[i];

    return x;
}

vector<double> testProgonka(const Matrix& slay, const vector<double>& x) {
    int n = slay.size();
    vector<double> e(n);
    for (int j = 0; j < n; j++) {
        double y = 0;
        for (int i = 0; i < n; i++) y += slay[j][i] * x[i];
        e[j] = round((slay[j][n] - y) * 1e6) / 1e6;
    }
    return e;
}

// ─── Норма матрицы / вектора (та же что в задаче 1) ─────────────────────────

double normaMatrix(const Matrix& B) {
    int n = B.size();
    double b1 = 0;
    for (int i = 0; i < n; i++) {
        double s = 0;
        for (int j = 1; j < (int)B[i].size() - 1; j++) s += fabs(B[i][j]);
        b1 = max(b1, s);
    }
    double b2 = 0;
    for (int j = 1; j < (int)B[0].size() - 1; j++) {
        double s = 0;
        for (int i = 0; i < n; i++) s += fabs(B[i][j]);
        b2 = max(b2, s);
    }
    double b3 = 0;
    for (int i = 0; i < n; i++)
        for (int j = 1; j < (int)B[i].size() - 1; j++)
            b3 += B[i][j] * B[i][j];
    b3 = sqrt(b3);
    return min(b1, min(b2, b3));
}

double normaVector(const vector<double>& v) {
    double b1 = 0, b2 = 0, b3 = 0;
    for (double val : v) {
        b1 = max(b1, fabs(val));
        b2 += fabs(val);
        b3 += val * val;
    }
    return min(b1, min(b2, sqrt(b3)));
}

// ─── Метод Зейделя ───────────────────────────────────────────────────────────

struct IterRow { int k; vector<double> x; double dx, ei; };

tuple<Matrix, vector<double>, vector<IterRow>> zeidelMethod(Matrix slay, double e) {
    int n = slay.size();

    for (int i = 0; i < n; i++) {
        slay[i].insert(slay[i].begin(), 0.0);
        swap(slay[i][i + 1], slay[i][0]);
    }
    for (int j = 0; j < n; j++)
        for (int i = 1; i < (int)slay[j].size() - 1; i++)
            slay[j][i] = -slay[j][i];
    for (int i = 0; i < n; i++)
        slay[i] = division(slay[i], slay[i][0]);

    vector<double> x0(n);
    for (int i = 0; i < n; i++) x0[i] = slay[i].back();

    double B = normaMatrix(slay);
    double ei = (B / (1.0 - B)) * normaVector(x0);

    vector<IterRow> iters;
    int k = 0;

    while (ei > e) {
        vector<double> x(n);
        for (int i = 0; i < n; i++) {
            double xi = 0;
            for (int j = 1; j < (int)slay[i].size() - 1; j++) {
                if (j - 1 < i)
                    xi += slay[i][j] * x[j - 1];
                else
                    xi += slay[i][j] * x0[j - 1];
            }
            xi += slay[i].back();
            x[i] = xi;
        }
        vector<double> diff(n);
        for (int i = 0; i < n; i++) diff[i] = x[i] - x0[i];
        double d = normaVector(diff);
        ei = (B / (1.0 - B)) * d;
        iters.push_back({k, x, d, ei});
        x0 = x;
        k++;
    }

    return {slay, x0, iters};
}

// ─── main ────────────────────────────────────────────────────────────────────

int main() {
    SetConsoleOutputCP(65001);   // UTF-8
    SetConsoleCP(65001);
    Matrix slay = {
        {8,  4,   0,  0,  11},
        {-2,  6,  -2,  0,   5},
        {0,  2, 9, -5,  15.5},
        {0,  0,  -2, 6,   9.5}
    };
    double e = 0.01;
    int n = slay.size();

    cout << fixed << setprecision(6);

    // ── Метод прогонки ──
    cout << string(50, '-') << "\nМетод прогонки\n" << string(50, '-') << "\n";
    auto xProg = progonka(slay);
    auto eProg = testProgonka(slay, xProg);

    cout << "Решение методом прогонки:\n";
    for (int i = 0; i < n; i++)
        cout << "x" << i + 1 << ": " << xProg[i] << "\n";
    cout << string(50, '-') << "\nПроверка:\n";
    for (int i = 0; i < n; i++)
        cout << "e" << i + 1 << ": " << eProg[i] << "\n";
    cout << string(50, '-') << "\n";

    // ── Метод Зейделя ──
    auto [newSlay, x0, iters] = zeidelMethod(slay, e);

    cout << "\n" << string(50, '-') << "\nМетод Зейделя\n" << string(50, '-') << "\n";
    cout << "Преобразованная матрица\n" << string(50, '-') << "\n";
    for (auto& row : newSlay) {
        for (double v : row) cout << setw(12) << v;
        cout << "\n";
    }

    cout << string(50, '-') << "\nРешение\n";
    cout << setw(4) << "k";
    for (int i = 0; i < n; i++) cout << setw(11) << ("x" + to_string(i + 1));
    cout << setw(11) << "||dx||" << setw(11) << "ei" << "\n";

    for (auto& r : iters) {
        cout << setw(4) << r.k;
        for (double v : r.x) cout << setw(11) << v;
        cout << setw(11) << r.dx << setw(11) << r.ei << "\n";
    }

    cout << string(50, '-') << "\nОтвет\n";
    auto& last = iters.back();
    cout << setw(4) << last.k;
    for (double v : last.x) cout << setw(11) << v;
    cout << setw(11) << last.dx << setw(11) << last.ei << "\n";

    return 0;
}
