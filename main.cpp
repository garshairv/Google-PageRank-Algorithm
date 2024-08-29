//Name: Misam Ibrahimi, Student#: A01204003
//Name: Amirgarsha Iravani, Student#: A01231321

#include <iostream>
#include "Matrix.hpp"
using namespace std;

int main() {
    vector<double> values = fileReader();
    Matrix *connectivityMatrix = new Matrix(values);
    Matrix *Q = new Matrix(values);
    generateImportanceMatrix(*Q);
    generateQMatrix(*connectivityMatrix);
    generateTransitionMatrix(*connectivityMatrix, *Q);
    markovProcess(*Q);
    delete(Q);
    delete(connectivityMatrix);
    return 0;
}
