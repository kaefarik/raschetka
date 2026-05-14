#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <string>
#include <windows.h>

using namespace std;

typedef vector<vector<double>> Matrix;

pair<vector<double>, Matrix> jacobiEigen(Matrix A, double e) {
    int n = A.size();

    // D — рабочая матрица (копия A), V — матрица собственных векторов (единичная)
    Matrix D = A;
    Matrix V(n, vector<double>(n, 0.0));
    for (int i = 0; i < n; i++) V[i][i] = 1.0;

    while (true) {
        // Находим максимальный внедиагональный элемент
        double maxVal = 0.0;
        int p = 0, q = 0;
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++)
                if (fabs(D[i][j]) > maxVal) {
                    maxVal = fabs(D[i][j]);
                    p = i; q = j;
                }

        if (maxVal < e) break;

        // Угол поворота
        double theta;
        if (D[p][p] == D[q][q])
            theta = M_PI / 4.0;
        else
            theta = 0.5 * atan(2.0 * D[p][q] / (D[p][p] - D[q][q]));

        double ct = cos(theta), st = sin(theta);

        // Обновляем D
        Matrix Dnew = D;
        Dnew[p][p] = D[p][p]*ct*ct + 2*D[p][q]*st*ct + D[q][q]*st*st;
        Dnew[q][q] = D[p][p]*st*st - 2*D[p][q]*st*ct + D[q][q]*ct*ct;
        Dnew[p][q] = Dnew[q][p] = 0.0;

        for (int i = 0; i < n; i++) {
            if (i != p && i != q) {
                Dnew[i][p] = Dnew[p][i] = D[i][p]*ct + D[i][q]*st;
                Dnew[i][q] = Dnew[q][i] = -D[i][p]*st + D[i][q]*ct;
            }
        }
        D = Dnew;

        // Обновляем V
        for (int i = 0; i < n; i++) {
            double vip = V[i][p], viq = V[i][q];
            V[i][p] =  vip*ct + viq*st;
            V[i][q] = -vip*st + viq*ct;
        }
    }

    // Собственные значения — диагональ D
    vector<double> eigenvalues(n);
    for (int i = 0; i < n; i++) eigenvalues[i] = D[i][i];

    // Собственные векторы — столбцы V (транспонируем для удобства вывода)
    Matrix eigenvectors(n, vector<double>(n));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            eigenvectors[i][j] = V[j][i];

    return {eigenvalues, eigenvectors};
}

int main() {
    SetConsoleOutputCP(65001);   // UTF-8
    SetConsoleCP(65001);
    Matrix A = {
        {68, 4,  3},
        { 4, 32, 4},
        { 3,  4, 84}
    };
    double e = 0.01;

    auto [eigenvalues, eigenvectors] = jacobiEigen(A, e);

    cout << fixed << setprecision(6);
    cout << "Собственные значения:\n";
    for (int i = 0; i < (int)eigenvalues.size(); i++)
        cout << "λ" << i + 1 << " = " << eigenvalues[i] << "\n";

    cout << string(50, '-') << "\n";
    cout << "\nСобственные векторы:\n";
    for (int i = 0; i < (int)eigenvectors.size(); i++) {
        cout << "v" << i + 1 << " = [";
        for (int j = 0; j < (int)eigenvectors[i].size(); j++) {
            cout << eigenvectors[i][j];
            if (j + 1 < (int)eigenvectors[i].size()) cout << ", ";
        }
        cout << "]\n";
    }

    return 0;
}
