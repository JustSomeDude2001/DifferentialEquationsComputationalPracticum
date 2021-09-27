#include <bits/stdc++.h>

using namespace std;

void euler(float (*y)(float), float (*f)(float, float), float y_0, float x_0, float h,
           //Output segment pointers:
           float (*x_res), float (*y_res), float (*y_euler), float (*lte), float (*gte),
           size_t n) {
    // Initial value assignment
    (*x_res) = x_0;
    (*y_res) = y_0;
    (*y_euler) = y_0;
    (*lte) = 0.0f;
    (*gte) = 0.0f;

    // The method implementation
    for (int i = 0; i < n - 1; i++) {
        x_res[i + 1] = x_res[i] + h;
        y_res[i + 1] = y(x_res[i + 1]);
        y_euler[i + 1] = y_euler[i] + h * f(x_res[i], y_euler[i]);
        lte[i + 1] = (y_res[i] + h * f(x_res[i], y_res[i])) -
                     (y_res[i + 1]);
        gte[i + 1] = y_euler[i + 1] - y_res[i + 1];
    }
}

void eulerImproved(float (*y)(float), float (*f)(float, float), float y_0, float x_0, float h,
           //Output segment pointers:
           float (*x_res), float (*y_res), float (*y_euler), float (*lte), float (*gte),
           size_t n) {
    // Initial value assignment
    (*x_res) = x_0;
    (*y_res) = y_0;
    (*y_euler) = y_0;
    (*lte) = 0.0f;
    (*gte) = 0.0f;

    // The method implementation
    for (int i = 0; i < n - 1; i++) {
        x_res[i + 1] = x_res[i] + h;
        y_res[i + 1] = y(x_res[i + 1]);
        y_euler[i + 1] = y_euler[i] + h * f(x_res[i] + h * 0.5f, y_euler[i] + h * 0.5f * f(x_res[i], y_euler[i]));
        lte[i + 1] = (y_res[i] + h * f(x_res[i] + h * 2.0f, y_res[i] + h * 2.0f * f(x_res[i], y_res[i]))) -
                     (y_res[i + 1]);
        gte[i + 1] = y_euler[i + 1] - y_res[i + 1];
    }
}

void rungeKutter(float (*y)(float), float (*f)(float, float), float y_0, float x_0, float h,
           //Output segment pointers:
           float (*x_res), float (*y_res), float (*y_euler), float (*lte), float (*gte),
           size_t n) {
    // Initial value assignment
    (*x_res) = x_0;
    (*y_res) = y_0;
    (*y_euler) = y_0;
    (*lte) = 0.0f;
    (*gte) = 0.0f;

    // The method implementation
    for (int i = 0; i < n - 1; i++) {
        x_res[i + 1] = x_res[i] + h;
        y_res[i + 1] = y(x_res[i + 1]);

        float k1 = f(x_res[i], y_euler[i]);
        float k2 = f(x_res[i] + h * 0.5f, y_euler[i] + h * k1 * 0.5f);
        float k3 = f(x_res[i] + h * 0.5f, y_euler[i] + h * k2 * 0.5f);
        float k4 = f(x_res[i] + h, y_euler[i] + h * k3);

        y_euler[i + 1] = y_euler[i] + h * (1.0f / 6.0f) * (k1 + 2.0f * k2 + 2.0f * k3 + k4);
        lte[i + 1] = (y_res[i] + h * f(x_res[i] + h * 2.0f, y_res[i] + h * 2.0f * f(x_res[i], y_res[i]))) -
                     (y_res[i + 1]);
        gte[i + 1] = y_euler[i + 1] - y_res[i + 1];
    }
}

/**void applyIterativeMethod(void (*method)(),
                          float (*y)(float), float (*f)(float, float), float y_0, float x_0, float h,
                          //Output segment pointers:
                          float (*x_res), float (*y_res), float (*y_euler), float (*lte), float (*gte),
                          size_t n) {
    
    // Initial value assignment
    (*x_res) = x_0;
    (*y_res) = y_0;
    (*y_euler) = y_0;
    (*lte) = 0.0f;
    (*gte) = 0.0f;

    // The method implementation
    for (int i = 0; i < n - 1; i++) {
        x_res[i + 1] = x_res[i] + h;
        y_res[i + 1] = y(x_res[i + 1]);
        y_euler[i + 1] = y_euler[i] + h * f(x_res[i], y_euler[i]);
        lte[i + 1] = (y_res[i] + h * f(x_res[i], y_res[i])) -
                     (y_res[i + 1]);
        gte[i + 1] = y_euler[i + 1] - y_res[i + 1];
    }
    
    
}**/

string columnNames[5] = {
    "X_result",
    "Y_result",
    "Y_method",
    "LTE",
    "GTE"
};

void outputDiffEqTable(float (*x_res), float (*y_res), float (*y_euler), float (*lte), float (*gte),
                       size_t n, int printDist, int cellWidth) {
    
    string div = "";

    // Adapting divisor
    for (int i = 0; i < printDist; i++) {
        div += ' ';
    }

    // Adapting column names and printing
    for (int i = 0; i < 5; i++) {
        while (columnNames[i].size() < cellWidth) {
            columnNames[i] += ' ';
        }
        cout << columnNames[i] << div;
    }
    
    cout << '\n';

    cout.precision(cellWidth);
    cout.width(cellWidth);
    
    // Printing objects
    for (int i = 0; i < n; i++) {
        cout << *(x_res + i) << div;
        cout << *(y_res + i) << div;
        cout << *(y_euler + i) << div;
        cout << *(lte + i) << div;
        cout << *(gte + i) << div;
    
        cout << '\n';
    }
}

float y_sharp(float x, float y) {
    return (y * y + x * y - x * x) / (x * x);
}

float y_exact(float x) {
    return (x * (1.0f + (x * x) / 3.0f)) / (1.0f - (x * x) / 3.0f);
}

int main() {
    int n;
    cin >> n;

    float h;
    cin >> h;

    float x_res[n * 100];
    float y_res[n * 100];
    float y_euler[n * 100];
    float lte[n * 100];
    float gte[n * 100];

    float ratio[n * 100];

    euler(&y_exact, &y_sharp, 2.0f, 1.0f, h, x_res, y_res, y_euler, lte, gte, n);

    outputDiffEqTable(x_res, y_res, y_euler, lte, gte, n, 5, 10);

    for (int i = 0; i < n; i++) {
        ratio[i] = lte[i];
    }

    euler(&y_exact, &y_sharp, 2.0f, 1.0f, h / 10.0f, x_res, y_res, y_euler, lte, gte, n * 10);

    outputDiffEqTable(x_res, y_res, y_euler, lte, gte, n, 5, 10);

    for (int i = 0; i < n; i++) {
        ratio[i] = ratio[i] / lte[i * 10];
        cout << ratio[i] << '\n';
    }

    /**
    eulerImproved(&y_exact, &y_sharp, 2.0f, 1.0f, h, x_res, y_res, y_euler, lte, gte, n);

    outputDiffEqTable(x_res, y_res, y_euler, lte, gte, n, 5, 10);
    
    rungeKutter(&y_exact, &y_sharp, 2.0f, 1.0f, h, x_res, y_res, y_euler, lte, gte, n);

    outputDiffEqTable(x_res, y_res, y_euler, lte, gte, n, 5, 10);
    **/
}