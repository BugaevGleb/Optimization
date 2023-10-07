#include <iostream>
#include <vector>
#include <float.h>
#include <cmath>
#include <iomanip>

using namespace std;

const double epsilon = 1e-10;

class Matrix {
private:
    int n;
    int m;
    double **matrix;

public:

    Matrix(int n, int m) {
        this->n = n;
        this->m = m;
        // Allocate memory for the matrix
        matrix = (double **) malloc(n * sizeof(double *));
        for (int i = 0; i < n; i++) {
            matrix[i] = (double *) malloc(m * sizeof(double));
        }
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                matrix[i][j] = 0.0;
            }
        }
    }

    // Copy constructor
    Matrix(Matrix &a) : Matrix(a.getN(), a.getM()) {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                matrix[i][j] = a[i][j];
            }
        }
    }

    // Destructor
    ~Matrix() {
        for (int i = 0; i < n; ++i) {
            delete matrix[i];
        }
        delete matrix;
    }

    int getN() { return n; }

    int getM() { return m; }

    double *operator[](int index) {
        return matrix[index];
    }

    vector<double> getColumn(int column) {
        vector<double> res(getN());
        for (int i = 0; i < getN(); ++i) {
            for (int j = 0; j < getM(); ++j) {
                if (j == column) {
                    res[i] = matrix[i][j];
                }
            }
        }
        return res;
    }


    friend istream &operator>>(istream &stream, Matrix &matrix1) {
        for (int i = 0; i < matrix1.getN(); ++i) {
            for (int j = 0; j < matrix1.getM(); ++j) {
                stream >> matrix1[i][j];
            }
        }
        return stream;

    }

    friend ostream &operator<<(ostream &stream, Matrix &matrix1) {

        for (int i = 0; i < matrix1.n; ++i) {
            for (int j = 0; j < matrix1.m; ++j) {
                if (fabs(matrix1[i][j]) < epsilon) matrix1[i][j] = 0.0;
                stream << matrix1[i][j] << " ";
            }
            stream << endl;
        }
        return stream;
    }
};

template<typename T>
void readVector(vector<T> &v, int n) {
    for (int i = 0; i < n; ++i) {
        T t;
        cin >> t;
        v.push_back(t);
    }
}


bool checkStandardFormVectorB(vector<double> &b) {
    for (double &i: b) {
        if (i < 0) {
            return false;
        }
    }
    return true;
}


int getIndexOfBasicVariable(Matrix &matrix, int column) {
    int pivotIndex = -1;
    for (int i = 0; i < matrix.getN(); ++i) {
        if (matrix[i][column] != 1 && matrix[i][column] != 0) return -1;
        if (matrix[i][column] == 1) {
            if (pivotIndex == -1) pivotIndex = i;
            else return -1;
        }
    }
    return pivotIndex;
}


double dotProduct(vector<double> &a, vector<double> &b) {
    if (a.size() != b.size()) return -1;
    double res = 0;
    for (int i = 0; i < a.size(); ++i) {
        res = res + a[i] * b[i];
    }
    return res;
}


// Works only for normal form. All constraints are less or equal and right side is non-negative
void solve(vector<double> &c, Matrix &matrix, vector<double> &b, bool isMinimize) {
    vector<double> c_b(matrix.getN(), 0);
    vector<int> basic;
    vector<double> z(matrix.getM());
    vector<double> fractions(matrix.getN());


    //1. Choose basic variables

    for (int column = 0; column < matrix.getM(); ++column) {
        int indexOfBasicVariable = getIndexOfBasicVariable(matrix, column);
        if (indexOfBasicVariable != -1) {
            basic.push_back(column);
            c_b[indexOfBasicVariable] = c[column];
        }
    }
    if (basic.size() != matrix.getN()) {
        cout << "The method is not applicable!\n";
        exit(0);
    }


    while (true) {
        int leadingColumn = -1;
        double maxValue = 0;
        //2. Calculate difference z_i and choose maximum positive z_i among them
        for (int i = 0; i < matrix.getM(); ++i) {
            vector<double> matrixColumn = matrix.getColumn(i);
            z[i] = dotProduct(c_b, matrixColumn) - c[i];
            if (z[i] > 0 && z[i] > maxValue) {
                maxValue = z[i];
                leadingColumn = i;
            }
        }
        if (leadingColumn == -1) {
            vector<double> res(b.size());
            double functionValue = dotProduct(c_b, b);
            cout << endl;
            for (int i = 0; i < b.size(); ++i) {
                if (basic[i] <= b.size() - 1) {
                    res[basic[i]] = b[i];
                }
            }
            for (int i = 0; i < b.size(); ++i) {
                cout << "x" << i + 1 << " = " << res[i] << endl;
            }
            cout << "Result of the function = " << functionValue * (isMinimize ? 1 : -1);
            break;
        }

        //3. Calculate fractions and choose minimum positive among them
        int leadingRow = -1;
        double minValue = DBL_MAX;
        for (int i = 0; i < matrix.getN(); ++i) {
            fractions[i] = b[i] / matrix[i][leadingColumn];
            if (fractions[i] > 0 && fractions[i] < minValue) {
                minValue = fractions[i];
                leadingRow = i;
            }
        }

        double leadingElement = matrix[leadingRow][leadingColumn];

        //4. Change one of basic variables and divide its row by leadingElement
        for (int i = 0; i < matrix.getM(); ++i) {

            matrix[leadingRow][i] /= leadingElement;
        }
        b[leadingRow] /= leadingElement;

        for (int i = 0; i < matrix.getN(); ++i) {

            for (int j = 0; j < matrix.getM(); ++j) {
                if (i != leadingRow && j != leadingColumn) {
                    matrix[i][j] = matrix[i][j] - matrix[i][leadingColumn] * matrix[leadingRow][j];
                }
            }
            if (i != leadingRow) b[i] = b[i] - matrix[i][leadingColumn] * b[leadingRow];
        }
        for (int i = 0; i < matrix.getN(); ++i) {
            if (i == leadingRow) {
                matrix[i][leadingColumn] = 1;
            } else {
                matrix[i][leadingColumn] = 0;
            }
        }

        // 4.2 Recalculate c_b
        c_b[leadingRow] = c[leadingColumn];
        basic[leadingRow] = leadingColumn;
    }
}


vector<double> checkLinearityAndConvertVector(vector<string> v) {
    vector<double> f;
    for (int i = 0; i < size(v); i++) {
        try {
            f.push_back(stod(v[i]));
        }
        catch (...) {
            cout << "The method is not applicable!\n";
            exit(0);
        }
    }
    return f;
}


int main() {
    cout << "Type 1 if you want to minimize function and 2 if maximize:\n";
    bool isMinimize;

    int temp;
    cin >> temp;
    isMinimize = temp == 1;


    cout << "Enter number of non-slack variables:\n";
    int numberOfVariables;
    cin >> numberOfVariables;

    cout << "Enter number of constraint functions:\n";
    int numberOfConstraintFunctions;
    cin >> numberOfConstraintFunctions;

    cout << "Enter number of slack variables:\n";
    int numberOfSlackVariables;
    cin >> numberOfSlackVariables;

    cout << "Enter a vector of coefficients of objective function:\n";
    vector<string> initC;
    readVector(initC, numberOfVariables);
    cout << "Enter a matrix of coefficients of constraint function:\n";

    Matrix A(numberOfConstraintFunctions, numberOfVariables + numberOfSlackVariables);
    cin >> A;

    cout << "Enter a vector of right-hand side numbers:\n";
    vector<string> initB;
    readVector(initB, numberOfConstraintFunctions);
    vector<double> c, b;
    c = checkLinearityAndConvertVector(initC);
    b = checkLinearityAndConvertVector(initB);
    if (!checkStandardFormVectorB(b)) {
        cout << "The method is not applicable!\n";
        return 0;
    }
    if (!isMinimize) {
        for (int i = 0; i < numberOfVariables; ++i) {
            c[i] = -c[i];
        }
    }
    // Add zeros for slack variables
    for (int i = 0; i < numberOfSlackVariables; ++i) {
        c.push_back(0);
    }

    cout << "Enter the desired number of decimal places:\n";
    int decimalPlaces;
    cin >> decimalPlaces;
    cout << fixed;
    cout << setprecision(decimalPlaces);
    solve(c, A, b, isMinimize);

    return 0;
}
