const double epsilon = 1e-10;

#include <iostream>
#include <vector>
#include <iomanip>
#include <float.h>
#include "cmath"

using namespace std;

#include "Matrix.cpp"
#include "SimplexMethod.cpp"
#include "InteriorMethod.cpp"


int main() {

    cout << "Type 1 if you want to minimize function and 2 if maximize:\n";
    bool isMinimize;

    int temp;
    cin >> temp;
    isMinimize = temp == 1;
    int numberOfVariables;
    int numberOfConstraintFunctions;


    cout << "Enter number of variables\n";
    cin >> numberOfVariables;

    cout << "Enter number of constraint functions\n";
    cin >> numberOfConstraintFunctions;


    cout << "Enter a vector of coefficients of objective function ("  << numberOfVariables << " variables)"
         << ":\n";
    vector<double> c;
    readVector(c, numberOfVariables);
    cout << "Enter a matrix of coefficients of constraint function (" << numberOfConstraintFunctions << "x"
         << numberOfVariables << "):\n";

    Matrix A(numberOfConstraintFunctions, numberOfVariables);
    cin >> A;

    cout << "Enter a vector of right-hand side numbers (" <<numberOfConstraintFunctions <<" numbers):\n";
    vector<double> b;
    readVector(b, numberOfConstraintFunctions);


    if (isMinimize) {
        for (int i = 0; i < numberOfVariables; ++i) {
            c[i] = -c[i];
        }
    }
    vector<double> x;
    cout << "Enter the initial solution (" << numberOfVariables << " variables)" << ":\n";
    readVector(x, numberOfVariables);

    cout << "Enter the desired number of decimal places:\n";
    int decimalPlaces;
    cin >> decimalPlaces;
    cout << fixed;
    cout << setprecision(decimalPlaces);
    cout << "alpha = 0.5" << endl;
    solve(c, A, b, x, 0.5, isMinimize);
    cout << "--------------------" << endl;
    cout << "alpha = 0.9" << endl;
    solve(c, A, b, x, 0.9, isMinimize);
    cout << "--------------------" << endl;
    cout << "Simplex";
    if (!checkStandardFormVectorB(b)) {
        cout << "The method is not applicable!\n";
        return 0;
    }
    solve_simplex(c, A, b, isMinimize);

    return 0;
}
