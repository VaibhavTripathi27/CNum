#include "main.cpp"
#include <iostream>
#include <cmath>


//numpy.dot(A,B, out, s1, s2, d1, d2, axis)
void dot(int* matrix1, int* matrix2, int* result, int* shape1, int* shape2, int ndims1, int ndims2, int axis) 
{
    // Check that the number of dimensions is valid
    if (ndims1 != 2 || ndims2 != 2) {
        std::cerr << "Error: only two-dimensional matrices are supported" << std::endl;
        return;
    }

    // Check that the matrices are compatible for multiplication (the number of columns of the first matrix must equal the number of rows of the second matrix)
    if (shape1[1] != shape2[0]) {
        std::cerr << "Error: matrices are not compatible for multiplication" << std::endl;
        return;
    }

    // Compute the dot product
    int m = shape1[0];
    int n = shape1[1];
    int p = shape2[1];
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < p; ++j) {
            int sum = 0;
            for (int k = 0; k < n; ++k) {
                sum += matrix1[i * n + k] * matrix2[k * p + j];
            }
            result[i * p + j] = sum;
        }
    }
}

/* int main() {
    // Define two matrices with shape 2 x 3 and 3 x 2
    int matrix1[] = {1, 2, 3, 4, 5, 6};
    int matrix2[] = {7, 8, 9, 10, 11, 12};
    int shape1[] = {2, 3};
    int shape2[] = {3, 2};

    // Create an array to store the result with shape 2 x 2
    int result[4];

    // Call the dot function
    dot(matrix1, matrix2, result, shape1, shape2, 2, 2, -1);

    // Print the result
    std::cout << "Result of matrix multiplication:" << std::endl;
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            std::cout << result[i * 2 + j] << ' ';
        }
        std::cout << std::endl;
    }

    return 0;
}
*/

//numpy.vdot(vector1, vector2, length)

double vdot(double* vector1, double* vector2, int length)
{
    // Check that the vectors are one-dimensional
    if (length <= 0) {
        std::cerr << "Error: vectors must be one-dimensional" << std::endl;
        return 0.0;
    }

    // Compute the dot product
    double result = 0.0;
    for (int i = 0; i < length; ++i) {
        result += vector1[i] * vector2[i];
    }
    return result;
}

/* int main() {
    // Define two vectors with length 3
    double vector1[] = {1.0, 2.0, 3.0};
    double vector2[] = {4.0, 5.0, 6.0};

    // Call the vdot function
    double result = vdot(vector1, vector2, 3);

    // Print the result
    std::cout << "Result of dot product: " << result << std::endl;

    return 0;
}
*/

//numpy.rank(tensor, shape, dimensions)

int rank(double* tensor, int* shape, int dimensions) {
    // Check that the tensor is at least one-dimensional
    if (dimensions <= 0) {
        std::cerr << "Error: tensor must have at least one dimension" << std::endl;
        return -1;
    }

    // Return the number of dimensions
    return dimensions;
}

/* int main() {
    // Define a tensor with shape (3, 4)
    double tensor[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0};
    int shape[] = {3, 4};

    // Call the rank function
    int result = rank(tensor, shape, 2);

    // Print the result
    std::cout << "Rank of tensor: " << result << std::endl;

    return 0;
}
*/

//numpy.trace(tensor, shape, dimensions)

double trace(double* tensor, int* shape, int dimensions)
{
    // Check that the tensor is two-dimensional
    if (dimensions != 2) {
        std::cerr << "Error: tensor must be two-dimensional" << std::endl;
        return -1;
    }

    // Check that the matrix is square
    if (shape[0] != shape[1]) {
        std::cerr << "Error: matrix must be square" << std::endl;
        return -1;
    }

    // Get the size of the matrix
    int size = shape[0];

    // Compute the trace of the matrix
    double trace = 0.0;
    for (int i = 0; i < size; i++) {
        trace += tensor[i * size + i];
    }

    return trace;
}

/* int main() {
    // Define a matrix with shape (3, 3)
    double tensor[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
    int shape[] = {3, 3};

    // Call the trace function
    double result = trace(tensor, shape, 2);

    // Print the result
    std::cout << "Trace of matrix: " << result << std::endl;

    return 0;
}
*/

//numpy.solve(A, B, N)

void solve(double* A, double* b, int N) {
    // Forward elimination
    for (int k = 0; k < N - 1; k++) {
        for (int i = k + 1; i < N; i++) {
            double factor = A[i * N + k] / A[k * N + k];
            b[i] -= factor * b[k];
            for (int j = k + 1; j < N; j++) {
                A[i * N + j] -= factor * A[k * N + j];
            }
        }
    }

    // Backward substitution
    for (int i = N - 1; i >= 0; i--) {
        for (int j = i + 1; j < N; j++) {
            b[i] -= A[i * N + j] * b[j];
        }
        b[i] /= A[i * N + i];
    }
}


/* int main() {
    // Define a linear system of equations
    double A[] = {2.0, 1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 2.0};
    double b[] = {1.0, 2.0, 3.0};
    int N = 3;

    // Call the solve function
    solve(A, b, N);

    // Print the solution
    std::cout << "Solution of linear system of equations:" << std::endl;
    for (int i = 0; i < N; i++) {
        std::cout << "x" << i << " = " << b[i] << std::endl;
    }

    return 0;
}
*/

//numpy.inv(A, n)

void inv(double* A, int N) {
    // Augment the identity matrix to A
    double* B = new double[N * 2 * N];
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            B[i * 2 * N + j] = A[i * N + j];
            B[i * 2 * N + j + N] = (i == j) ? 1.0 : 0.0;
        }
    }

    // Forward elimination
    for (int k = 0; k < N - 1; k++) {
        for (int i = k + 1; i < N; i++) {
            double factor = B[i * 2 * N + k] / B[k * 2 * N + k];
            for (int j = k + 1; j < 2 * N; j++) {
                B[i * 2 * N + j] -= factor * B[k * 2 * N + j];
            }
        }
    }

    // Backward substitution
    for (int i = N - 1; i >= 0; i--) {
        for (int j = i + 1; j < N; j++) {
            for (int k = N; k < 2 * N; k++) {
                B[i * 2 * N + k] -= B[i * 2 * N + j] * B[j * 2 * N + k];
            }
        }
        for (int j = N; j < 2 * N; j++) {
            B[i * 2 * N + j] /= B[i * 2 * N + i];
        }
    }

    // Copy the inverse matrix to A
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            A[i * N + j] = B[i * 2 * N + j + N];
        }
    }

    delete[] B;
}

/* int main() {
    // Define a square matrix
    double A[] = {2.0, 1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 2.0};
    int N = 3;

    // Call the inv function
    inv(A, N);

    // Print the inverse matrix
    std::cout << "Inverse of square matrix:" << std::endl;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            std::cout << A[i * N + j] << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}
*/



double determinant(double** arr, int n) {
    double det = 0;
    if (n == 2) {
        det = arr[0][0] * arr[1][1] - arr[0][1] * arr[1][0];
        return det;
    }
    double** temp = new double*[n - 1];
    for (int i = 0; i < n - 1; i++) {
        temp[i] = new double[n - 1];
    }
    for (int i = 0; i < n; i++) {
        int subindex = 0;
        for (int j = 1; j < n; j++) {
            int subindex2 = 0;
            for (int k = 0; k < n; k++) {
                if (k == i) {
                    continue;
                }
                temp[subindex][subindex2] = arr[j][k];
                subindex2++;
            }
            subindex++;
        }
        det += pow(-1, i) * arr[0][i] * determinant(temp, n - 1);
    }
    for (int i = 0; i < n - 1; i++) {
        delete[] temp[i];
    }
    delete[] temp;
    return det;
}

/* int main() {
    int n;
    std::cout << "Enter the dimension of the matrix: ";
    std::cin >> n;
    double** arr = new double*[n];
    for (int i = 0; i < n; i++) {
        arr[i] = new double[n];
    }
    std::cout << "Enter the elements of the matrix: " << std::endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            std::cin >> arr[i][j];
        }
    }
    double det = determinant(arr, n);
    std::cout << "The determinant of the matrix is: " << det << std::endl;
    for (int i = 0; i < n; i++) {
        delete[] arr[i];
    }
    delete[] arr;
    return 0;
}
*/


