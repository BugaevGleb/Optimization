

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

void solve_simplex(vector<double> &c, Matrix &matrix, vector<double> &b, bool isMinimize) {
    vector<double> c_b(matrix.getN(), 0);
    vector<int> basic(matrix.getN());
    int numberOfBasicVariables = 0;
    vector<double> z(matrix.getM());
    vector<double> fractions(matrix.getN());


    //1. Choose basic variables

    for (int column = 0; column < matrix.getM(); ++column) {
        int indexOfBasicVariable = getIndexOfBasicVariable(matrix, column);
        if (indexOfBasicVariable != -1) {
            numberOfBasicVariables++;
            basic[indexOfBasicVariable] = column;
            c_b[indexOfBasicVariable] = c[column];
        }
    }
    if (numberOfBasicVariables != matrix.getN()) {
        cout << endl << "The method is not applicable!\n";
        exit(0);
    }


    while (true) {
        int leadingColumn = -1;
        double maxValue = 0;
        //2. Calculate difference z_i and choose maximum positive z_i among them
        for (int i = 0; i < matrix.getM(); ++i) {
            vector<double> matrixColumn = matrix.getColumn(i);
            z[i] = dotProduct(c_b, matrixColumn) - c[i];
            if (z[i] < 0 && z[i] < maxValue) {
                maxValue = z[i];
                leadingColumn = i;
            }
        }
        if (leadingColumn == -1) {
            vector<double> res(matrix.getM());
            double functionValue = dotProduct(c_b, b);
            cout << endl;
            for (int i = 0; i < res.size(); ++i) {
                if (i < basic.size()) {
                    res[basic[i]] = b[i];
                }
            }
            for (int i = 0; i < res.size(); ++i) {
                cout << "x" << i + 1 << " = " << res[i] << endl;
            }
            cout << "Value of the function = " << functionValue * (isMinimize ? -1 : 1);
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
        if(leadingRow==-1){
            cout<<endl<<"The method is not applicable!"<<endl;
            exit(0);
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