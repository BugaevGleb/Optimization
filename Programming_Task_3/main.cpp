#include <iostream>
#include <vector>
#include <float.h>
#include "cmath"

using namespace std;

class Matrix {
protected:
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
            //  delete matrix[i];
        }
        // delete matrix;
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

    // Return the index of row (below current row) with maximum absolute value in current column

    int findRowNumberWithMaximumAbsolutValueInColumnBelow(int column) {
        double maxValue = matrix[column][column];
        int resRow = column;
        for (int row = column + 1; row < n; ++row) {
            if (abs(matrix[row][column]) > abs(maxValue)) {
                maxValue = abs(matrix[row][column]);
                resRow = row;
            }
        }
        return resRow;
    }

    double getMinimalNegativeValue() {
        double res = DBL_MIN;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                if (matrix[i][j] < 0)
                    res = min(res, matrix[i][j]);
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
                if (fabs(matrix1[i][j]) < 0.00001) matrix1[i][j] = 0.0;
                stream << matrix1[i][j] << " ";
            }
            stream << endl;
        }
        return stream;
    }

    Matrix operator-(Matrix &input) {
        Matrix resMatrix(n, m);
        if (this->n == input.n && this->m == input.m) {

            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < m; ++j) {
                    resMatrix[i][j] = this->matrix[i][j] - input[i][j];
                }
            }
        } else {
            cout << "Error: the dimensional problem occurred. Please check the dimensions of the input" << endl;
            exit(1);

        }
        return resMatrix;
    }

    Matrix operator+(Matrix &input) {
        Matrix resMatrix(n, m);
        if (this->n == input.n && this->m == input.m) {
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < m; ++j) {
                    resMatrix[i][j] = input[i][j] + this->matrix[i][j];
                }
            }
        } else {
            cout << "Error: the dimensional problem occurred. Please check the dimensions of the input" << endl;
            exit(1);
        }
        return resMatrix;
    }

    Matrix &operator*(Matrix &secondMatrix) {
        auto resMatrix = new Matrix{n, secondMatrix.m};
        if (m == secondMatrix.n) {
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < secondMatrix.m; ++j) {
                    double res = 0.0;
                    for (int k = 0; k < m; ++k) {
                        res = res + matrix[i][k] * secondMatrix[k][j];

                    }
                    (*resMatrix)[i][j] = res;

                }
            }
        } else {
            cout << "Error: the dimensional problem occurred. Please check the dimensions of the input" << endl;
            exit(1);
        }
        return *resMatrix;
    }

    Matrix &operator*(double c) {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                matrix[i][j] = matrix[i][j] * c;

            }
        }
        return *this;
    }
};


class SquareMatrix : public Matrix {
public:
    SquareMatrix(int n) : Matrix(n, n) {}

    SquareMatrix(Matrix &a) : Matrix(a) {

    }

    int getSize() {
        return n;
    }

    using Matrix::operator=;
};

class IdentityMatrix : public SquareMatrix {
public:
    IdentityMatrix(int n) : SquareMatrix(n) {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i == j) matrix[i][j] = 1.0;
                else matrix[i][j] = 0.0;
            }
        }
    }

    using Matrix::operator=;

};

class EliminationMatrix : public IdentityMatrix {
public:
    EliminationMatrix(SquareMatrix &a, int i, int j) : IdentityMatrix(a.getSize()) {
        matrix[i][j] = -a[i][j] / a[j][j];

    }
};

class PermutationMatrix : public IdentityMatrix {
public:
    PermutationMatrix(int n, int row1, int row2) : IdentityMatrix(n) {

        swap(matrix[row1], matrix[row2]);

    }
};

class ColumnVector : public Matrix {
public:
    ColumnVector(int n) : Matrix(n, 1) {

    }

    ColumnVector(Matrix &a) : Matrix(a) {

    }

    ColumnVector(vector<double> a) : Matrix(a.size(), 1) {

        for (int i = 0; i < a.size(); ++i) {
            matrix[i][0] = a[i];

        }
    }

    // Returns the length of a ColumnVector
    double computeNorm() {
        double res = 0.0;
        for (int i = 0; i < n; ++i) {
            res += (*this)[i] + (*this)[i];
        }
        return ::sqrt(res);
    }

    // Overloaded operators
    double &operator[](int index) {
        return matrix[index][0];
    }

    using Matrix::operator=;
    using Matrix::operator*;
    using Matrix::operator+;

};

class DiagonalMatrix : public Matrix {
public:
    DiagonalMatrix(ColumnVector &a) : Matrix(a.getN(), a.getN()) {
        int n = a.getN();
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i == j) {
                    matrix[i][j] = a[i];
                } else {
                    matrix[i][j] = 0;
                }
            }
        }
    }

};


Matrix NorthWestCorner(ColumnVector S, Matrix C, ColumnVector D) {
    int x = 0, y = 0;
    Matrix supplied = Matrix(C.getN(), C.getM());
    while (y != C.getN() and x != C.getM()) {
        double temp = min(S[y], D[x]);
        S[y] -= temp;
        D[x] -= temp;
        supplied[y][x] = temp;
        if (S[y] == 0 and y < C.getN()) {
            y += 1;
        }
        if (D[x] == 0 and x < C.getM()) {
            x += 1;
        }
    }
    return supplied;
}


vector<double> analyzeRowAndColumnDiff(Matrix C, vector<bool> rowReady, vector<bool> columnReady) {
    auto row_diff = DBL_MIN, col_diff = DBL_MIN;
    int row_number = -1, col_number = -1;
    double temp1, temp2;
    for (int y = 0; y < C.getN(); y++) {
        temp1 = DBL_MAX;
        temp2 = DBL_MAX;
        if (!rowReady[y]) {
            continue;
        }
        for (int x = 0; x < C.getM(); x++) {
            if (!columnReady[x]) {
                continue;
            }
            if (C[y][x] < temp1) {
                temp2 = temp1;
                temp1 = C[y][x];
            } else if (C[y][x] < temp2) {
                temp2 = C[y][x];
            }
        }
        if (abs(temp2 - temp1) > row_diff) {
            row_diff = abs(temp2 - temp1);
            row_number = y;
        }
    }
    for (int x = 0; x < C.getM(); x++) {
        temp1 = DBL_MAX;
        temp2 = DBL_MAX;
        if (!columnReady[x]) {
            continue;
        }
        for (int y = 0; y < C.getN(); y++) {
            if (!rowReady[y]) {
                continue;
            }
            if (C[y][x] < temp1) {
                temp2 = temp1;
                temp1 = C[y][x];
            } else if (C[y][x] < temp2) {
                temp2 = C[y][x];
            }
        }
        if (temp2 - temp1 > col_diff) {
            col_diff = temp2 - temp1;
            col_number = x;
        }
    }
    vector<double> answer = {row_diff, col_diff, (double) row_number, (double) col_number};
    return answer;
}

Matrix VogelApproximation(ColumnVector S, Matrix C, ColumnVector D) {
    vector<bool> rowReady;
    vector<bool> columnReady;
    Matrix supplied = Matrix(C.getN(), C.getM());
    for (int i = 0; i < C.getN(); i++) {
        rowReady.push_back(true);
    }
    for (int i = 0; i < C.getM(); i++) {
        columnReady.push_back(true);
    }
    while (true) {
        int min_x, min_y;
        double row_diff, col_diff;
        int row_number, col_number;
        vector<double> stats = analyzeRowAndColumnDiff(C, rowReady, columnReady);
        row_diff = stats[0];
        col_diff = stats[1];
        row_number = (int) stats[2];
        col_number = (int) stats[3];
        if (row_number == -1 and col_number == -1) {
            return supplied;
        } else {
            if (row_diff >= col_diff) {
                min_y = row_number;
                min_x = -1;
                for (int x = 0; x < C.getM(); x++) {
                    if (min_x == -1) {
                        if (columnReady[x]) {
                            min_x = x;
                        }
                    } else if (C[min_y][x] < C[min_y][min_x] && columnReady[x]) {
                        min_x = x;
                    }
                }
            } else {
                min_x = col_number;
                min_y = -1;
                for (int y = 0; y < C.getN(); y++) {
                    if (min_y == -1) {
                        if (rowReady[y]) {
                            min_y = y;
                        }
                    } else if (C[y][min_x] < C[min_y][min_x] && rowReady[y]) {
                        min_y = y;
                    }
                }
            }
        }
        double temp = min(S[min_y], D[min_x]);
        S[min_y] -= temp;
        D[min_x] -= temp;
        supplied[min_y][min_x] = temp;
        if (S[min_y] == 0) {
            rowReady[min_y] = false;
        }
        if (D[min_x] == 0) {
            columnReady[min_x] = false;
        }
    }
}

Matrix RusselApproximation(ColumnVector S, Matrix C, ColumnVector D) {
    vector<bool> rowReady;
    vector<bool> columnReady;
    Matrix supplied = Matrix(C.getN(), C.getM());
    for (int i = 0; i < C.getN(); i++) {
        rowReady.push_back(true);
    }
    for (int i = 0; i < C.getM(); i++) {
        columnReady.push_back(true);
    }
    while (true) {
        vector<double> u;
        vector<double> v;
        for (int i = 0; i < C.getN(); i++) {
            u.push_back(-1);
        }
        for (int i = 0; i < C.getM(); i++) {
            v.push_back(-1);
        }
        for (int y = 0; y < C.getN(); y++) {
            for (int x = 0; x < C.getM(); x++) {
                if (!(rowReady[y] && columnReady[x])) {
                    continue;
                }
                if (u[y] < C[y][x]) {
                    u[y] = C[y][x];
                }
                if (v[x] < C[y][x]) {
                    v[x] = C[y][x];
                }
            }
        }
        Matrix deltas = Matrix(C);
        double max_neg = 0;
        vector<int> coords = {-1, -1};
        for (int y = 0; y < C.getN(); y++) {
            for (int x = 0; x < C.getM(); x++) {
                deltas[y][x] = C[y][x] - u[y] - v[x];
                if (deltas[y][x] < 0 && max_neg > deltas[y][x] && rowReady[y] && columnReady[x]) {
                    max_neg = deltas[y][x];
                    coords = {y, x};
                }
            }
        }
        if (max_neg == 0) {
            break;
        }
        int y = coords[0];
        int x = coords[1];
        double temp = min(S[y], D[x]);
        S[y] -= temp;
        D[x] -= temp;
        supplied[y][x] = temp;
        if (S[y] == 0) {
            rowReady[y] = false;
        }
        if (D[x] == 0) {
            columnReady[x] = false;
        }
    }
    return supplied;
}


void print_initial(ColumnVector S, Matrix C, ColumnVector D) {
    cout << "\t";
    for (int i = 0; i < D.getN(); i++) {
        cout << i + 1 << "\t";
    }
    cout << "supply\n";
    for (int i = 0; i < S.getN(); i++) {
        cout << i + 1 << "\t";
        for (int j = 0; j < C.getM(); j++) {
            cout << C[i][j] << "\t";
        }
        cout << S[i] << "\n";
    }
    cout << "demand\t";
    for (int i = 0; i < D.getN(); i++) {
        cout << D[i] << "\t";
    }
    cout << "\n\n";
}


double computeZ(Matrix solution, Matrix C) {
    double Z = 0;
    for (int i = 0; i < solution.getN(); i++) {
        for (int j = 0; j < solution.getM(); j++) {
            Z += solution[i][j] * C[i][j];
        }
    }
    return Z;
}


bool checkBalance(ColumnVector S, ColumnVector D) {
    double sum1 = 0, sum2 = 0;
    for (int i = 0; i < S.getN(); i++) {
        sum1 += S[i];
    }
    for (int i = 0; i < D.getN(); i++) {
        sum2 += D[i];
    }
    return sum1 == sum2;
}

bool checkApplicability(ColumnVector S, ColumnVector D, Matrix C) {
    for (int i = 0; i < S.getN(); i++) {
        if (S[i] < 0) {
            return false;
        }
    }
    for (int i = 0; i < D.getN(); i++) {
        if (D[i] < 0) {
            return false;
        }
    }
    for (int i = 0; i < C.getN(); i++) {
        for (int j = 0; j < C.getM(); j++) {
            if (C[i][j] < 0) {
                return false;
            }
        }
    }
    return true;
}


int main() {
    int n = 3, m = 4;
    ColumnVector S1(n);
    ColumnVector D1(m);
    Matrix C(n, m);
    cout << "Enter source vector:\n";
    cin >> S1;
    cout << "Enter cost matrix:\n";
    cin >> C;
    cout << "Enter demand vector:\n";
    cin >> D1;
    print_initial(S1, C, D1);
    if (!checkApplicability(S1, D1, C)) {
        cout << "This method is not applicable!";
        return 0;
    }
    if (!checkBalance(S1, D1)) {
        cout << "The problem is not balanced!";
        return 0;
    }
    Matrix answer1 = NorthWestCorner(S1, C, D1);
    Matrix answer2 = VogelApproximation(S1, C, D1);
    Matrix answer3 = RusselApproximation(S1, C, D1);
    cout << "North-West Corner rule:\n" << answer1 << "Z = " << computeZ(answer1, C) << "\n\n";
    cout << "Vogel's Approximation:\n" << answer2 << "Z = " << computeZ(answer2, C) << "\n\n";
    cout << "Russel's Approximation:\n" << answer3 << "Z = " << computeZ(answer3, C) << "\n\n";
}
