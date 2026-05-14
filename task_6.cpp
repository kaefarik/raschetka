#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <string>
#include <windows.h>

using namespace std;

struct SplineResult {
    vector<double> y, b, c, d;
};

SplineResult solveSpline(const vector<pair<double, double>>& pts) {
    int n = pts.size() - 1;

    vector<double> h(n), yv(pts.size());
    for (int i = 0; i < n; i++) {
        h[i]  = pts[i + 1].first - pts[i].first;
        yv[i] = pts[i].second;
    }
    yv[n] = pts[n].second;

    vector<double> alpha(n, 0.0);
    for (int i = 1; i < n; i++)
        alpha[i] = (3.0 / h[i])  * (yv[i + 1] - yv[i])
                 - (3.0 / h[i-1]) * (yv[i]     - yv[i-1]);

    vector<double> l(n + 1, 1.0), mu(n + 1, 0.0), z(n + 1, 0.0);

    for (int i = 1; i < n; i++) {
        l[i]  = 2.0 * (pts[i + 1].first - pts[i - 1].first) - h[i - 1] * mu[i - 1];
        mu[i] = h[i] / l[i];
        z[i]  = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
    }

    vector<double> b(n, 0.0), c(n + 1, 0.0), d(n, 0.0);

    for (int j = n - 1; j >= 0; j--) {
        c[j] = z[j] - mu[j] * c[j + 1];
        b[j] = (yv[j + 1] - yv[j]) / h[j] - h[j] * (c[j + 1] + 2.0 * c[j]) / 3.0;
        d[j] = (c[j + 1] - c[j]) / (3.0 * h[j]);
    }

    return {yv, b, c, d};
}

int main() {
    SetConsoleOutputCP(65001);   // UTF-8
    SetConsoleCP(65001);
    vector<pair<double, double>> points = {
        {1.415, 0.888551},
        {1.420, 0.889559},
        {1.425, 0.890637},
        {1.430, 0.891667},
        {1.435, 0.892687},
        {1.440, 0.893698}
    };

    auto [yv, b, c, d] = solveSpline(points);

    cout << fixed << setprecision(5);

    for (int i = 0; i < (int)points.size() - 1; i++) {
        double xi = points[i].first;
        cout << "Segment [" << xi << ", " << points[i + 1].first << "]:\n";
        cout << "S_" << i + 1 << "(x) = "
             << yv[i] << " + "
             << b[i]  << "(x - " << xi << ") + "
             << c[i]  << "(x - " << xi << ")^2 + "
             << d[i]  << "(x - " << xi << ")^3\n";
        cout << string(10, '-') << "\n";
    }

    return 0;
}
