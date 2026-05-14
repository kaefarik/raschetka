#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <string>
#include <windows.h>

using namespace std;

typedef vector<double> Row;
typedef vector<Row> Matrix;

// ─── Вспомогательные операции со строками ───────────────────────────────────

Row division(Row line, double n) {
    for (auto& v : line) v /= n;
    return line;
}

Row multiplication(Row line, double n) {
    for (auto& v : line) v *= n;
    return line;
}

Row addition(Row a, const Row& b) {
    for (int i = 0; i < (int)a.size(); i++) a[i] += b[i];
    return a;
}

// ─── Метод Гаусса ────────────────────────────────────────────────────────────

pair<Matrix, vector<double>> gauss(Matrix slay) {
    int n = slay.size();

    for (int i = 0; i < n; i++) {
        if (slay[i][i] == 0.0)
            swap(slay[i], slay[n - 1]);

        slay[i] = division(slay[i], slay[i][i]);

        for (int j = i + 1; j < n; j++)
            slay[j] = addition(slay[j], multiplication(slay[i], -slay[j][i]));
    }

    vector<double> x(n);
    x[n - 1] = slay[n - 1][n];

    for (int i = n - 2; i >= 0; i--) {
        x[i] = slay[i][n];
        for (int j = i + 1; j < n; j++)
            x[i] -= slay[i][j] * x[j];
    }

    return {slay, x};
}

vector<double> testGauss(const Matrix& orig, const vector<double>& x) {
    int n = orig.size();
    vector<double> e(n);
    for (int j = 0; j < n; j++) {
        double y = 0;
        for (int i = 0; i < n; i++) y += orig[j][i] * x[i];
        e[j] = round((orig[j][n] - y) * 1000.0) / 1000.0;
    }
    return e;
}

// ─── Метод простой итерации ──────────────────────────────────────────────────

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
    for (double x : v) {
        b1 = max(b1, fabs(x));
        b2 += fabs(x);
        b3 += x * x;
    }
    b3 = sqrt(b3);
    return min(b1, min(b2, b3));
}

struct IterRow { int k; double x1, x2, x3, dx, ei; };

tuple<Matrix, vector<double>, vector<IterRow>> simpleIteration(Matrix slay, double e) {
    int n = slay.size();

    // Преобразование: выносим диагональный элемент вперёд
    for (int i = 0; i < n; i++) {
        slay[i].insert(slay[i].begin(), 0.0);
        swap(slay[i][i + 1], slay[i][0]);
    }
    // Меняем знаки внедиагональных коэффициентов
    for (int j = 0; j < n; j++)
        for (int i = 1; i < (int)slay[j].size() - 1; i++)
            slay[j][i] = -slay[j][i];
    // Делим на диагональный
    for (int i = 0; i < n; i++)
        slay[i] = division(slay[i], slay[i][0]);

    vector<double> x0(n);
    for (int i = 0; i < n; i++) x0[i] = slay[i].back();

    double B = normaMatrix(slay);
    double ei = (B / (1.0 - B)) * normaVector(x0);

    vector<IterRow> iterations;
    int k = 0;

    while (ei > e) {
        vector<double> x(n);
        for (int i = 0; i < n; i++) {
            double xi = 0;
            for (int j = 1; j < (int)slay[i].size() - 1; j++)
                xi += slay[i][j] * x0[j - 1];
            xi += slay[i].back();
            x[i] = xi;
        }
        vector<double> diff(n);
        for (int i = 0; i < n; i++) diff[i] = x[i] - x0[i];

        double d = normaVector(diff);
        ei = (B / (1.0 - B)) * d;
        iterations.push_back({k, x[0], x[1], x[2], d, ei});
        x0 = x;
        k++;
    }

    return {slay, x0, iterations};
}

// ─── main ────────────────────────────────────────────────────────────────────

int main() {
    SetConsoleOutputCP(65001);   // UTF-8
    SetConsoleCP(65001);
    Matrix slay = {
        {1.61, 0.71, -0.05, -0.44},
        {-1.03,  -2.05, 0.87,  -1.16},
        { 2.50,  3.12, -6.03,  -7.50}
    };
    double e = 0.01;

    cout << fixed << setprecision(6);

    // ── Метод Гаусса ──
    auto [newSlay, x] = gauss(slay);
    auto eGauss = testGauss(slay, x);

    cout << string(50, '-') << "\nМетод Гаусса\n" << string(50, '-') << "\n";
    cout << "Преобразованная матрица\n" << string(50, '-') << "\n";
    for (auto& row : newSlay) {
        for (double v : row) cout << setw(12) << v;
        cout << "\n";
    }
    cout << string(50, '-') << "\nРешение\n";
    for (int i = 0; i < (int)x.size(); i++)
        cout << "x" << i + 1 << " : " << x[i] << "\n";
    cout << string(50, '-') << "\nПроверка для метода Гаусса\n" << string(50, '-') << "\n";
    for (int i = 0; i < (int)eGauss.size(); i++)
        cout << "e для x" << i + 1 << " = " << eGauss[i] << "\n";
    cout << string(50, '-') << "\n";

    // ── Метод простой итерации ──
    auto [newSlayIt, x0, iters] = simpleIteration(slay, e);

    cout << "\n" << string(50, '-') << "\nМетод простой итерации\n" << string(50, '-') << "\n";
    cout << "Преобразованная матрица\n" << string(50, '-') << "\n";
    for (auto& row : newSlayIt) {
        for (double v : row) cout << setw(12) << v;
        cout << "\n";
    }
    cout << string(50, '-') << "\nРешение\n";
    cout << setw(4) << "k"
         << setw(11) << "x1"
         << setw(11) << "x2"
         << setw(11) << "x3"
         << setw(11) << "||dx||"
         << setw(11) << "ei" << "\n";
    for (auto& r : iters)
        cout << setw(4) << r.k
             << setw(11) << r.x1
             << setw(11) << r.x2
             << setw(11) << r.x3
             << setw(11) << r.dx
             << setw(11) << r.ei << "\n";
    cout << string(50, '-') << "\n";
    for (int i = 0; i < (int)x0.size(); i++)
        cout << "x" << i + 1 << " = " << x0[i] << "\n";
    cout << string(50, '-') << "\n";

    return 0;
}
