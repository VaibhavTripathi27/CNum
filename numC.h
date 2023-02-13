//This library is disributed in the hope that it would be useful,
//under the pretext of being a library similar to numpy in python.
//This library is a free software, and can be redistributed and/or
//modified as Free Open Source Software. This is not an internal header file,
//and the source code is available at :
//{hyperlink}

//#include "functions.cpp"
//#include "linalg.cpp"
#include <array>
#include <iostream> 
#include <sstream>
#include <string>
#include <algorithm>
#include <cmath>
#include <climits>
#include <cstring>
#include <vector>

using namespace std;


template <typename dType,
          int dims> // Generic datatype dType and the dimensions of the tensor

class tensor {
private:
  int itemsize, nbytes, stride = 0, dimensions = dims, shape_i[dims], size_i;

  dType *data; // The data is stored in a 1-D flexible array member

public:
  int size, shape[dims], indexing[dims], ndim = dims;
  bool Sliced = false;

  // Constructors
  tensor(dType data_[], int shape_[]) {
    // Size of each element of the tensor
    itemsize = sizeof(dType);
    // Calculating the size and storing the shape
    size = 1;
    for (int i = 0; i < ndim; i++) {
      shape[i] = shape_[i];
      shape_i[i] = shape_[i];
      size *= shape_[i];
    }

    nbytes = itemsize * size; // Total size of the tensor
    size_i = size;
    // Making the indexing array for faster access
    indexing[ndim - 1] = 1;
    for (int i = ndim - 2, temp = 1; i >= 0; i--) {
      indexing[i] = temp * shape[i + 1];
      temp = indexing[i];
    }

    // Storing the data
    data = new dType[size];
    for (int i = 0; i < size; i++) {
      data[i] = data_[i];
    }
  }

  // Copy Constructor
  tensor(const tensor &obj) {
    itemsize = obj.itemsize, nbytes = obj.nbytes, stride = obj.stride,
    size_i = obj.size_i;
    dimensions = obj.dimensions;
    size = obj.size, ndim = obj.ndim;
    Sliced = obj.Sliced;

    for (int i = 0; i < dimensions; i++) {
      shape_i[i] = obj.shape_i[i];
      // cout<<"i "<<shape_i[i]<<endl;
      shape[i] = obj.shape[i];
      // cout<<"s "<<shape[i]<<endl;
      indexing[i] = obj.indexing[i];
      // cout<<"indx "<<indexing[i]<<endl;
    }

    // Copying the Data
    data = new dType[size_i];
    for (int i = 0; i < size_i; i++) {
      data[i] = obj.data[i];
    }
  }

  // Constructor with stride variable
  tensor(dType data_[], int shape_[], int pdim, int indexing_arr[],
         int stri) { // pdim = dims of parent tensor, indexing_arr = indexing
                     // array of parent tensor, stri = stride
    // Size of each element of the tensor
    itemsize = sizeof(dType);
    // Calculating the size and storing the shape
    size = 1;
    for (int i = 0; i < ndim; i++) {
      shape[i] = shape_[i];
      size *= shape_[i];
    }
    nbytes = itemsize * size; // Total size of the tensor

    // Making the indexing array for faster access
    for (int i = 0; i < pdim; i++) {
      indexing[i] = indexing_arr[pdim - ndim + i];
    }

    // Storing the data
    data = new dType[size];
    for (int i = 0; i < size; i++) {
      data[i] = data_[i];
    }
    stride += stri; // Setting up the stride
  }

  // Manipulation Routiens
  int slice(int slicing_array[dims][2]) { // returns the size of the slice
    int new_dim = 0;
    bool flag = true; // For checking if some dimesions are reduced
    string shape_s = "";

    // Figure out the dimensions and shape of the new tensor
    for (int i = 0; i < dims; i++) {
      char s = '-';
      int dif = slicing_array[i][1] - slicing_array[i][0];

      if ((dif == 0) && flag) {
        // If only 1 of ith dimension is taken and no multiple of any dim has
        // been taken then assume it as gone
      } else {
        flag = false; // Multiple of a dim taken
        shape_s += to_string((dif + 1)) + s;
        new_dim++;
      }
    }

    // Compute the size, stride and shape array
    int str = 0;
    int siz = 1;
    stringstream ss(shape_s);
    string s;
    int i = 0;
    while (getline(ss, s, '-')) {
      shape[i] = stoi(s);
      siz *= shape[i];
      i++;
    }
    for (int i = 0; i < dimensions; i++)
      str += indexing[i] * slicing_array[i][0];

    // After this the first ndim elements of the shape array will be modified
    // cout<<"stride "<<str<<endl;
    //  Setting the state vars
    stride = str;
    ndim = new_dim;
    Sliced = true;
    size = siz;

    return siz;
  }

  void resetSlice() {
    stride = 0;
    ndim = dimensions;
    Sliced = false;
    for (int i = 0; i < dimensions; i++) {
      shape[i] = shape_i[i];
    }
    size = size_i;
  }

  // Overloaded Operators
  dType operator[](int loc[]) {
    int index = stride;
    // varible stride
    for (int i = 0; i < ndim; i++) {
      index +=
          loc[i] * indexing[dimensions - ndim +
                            i]; // Taking the last ndim number of coefficients
    }
    // cout<<"\nIndex =  "<<index<<" "<<stride<<endl;
    return data[index];
  }

  auto operator=(tensor &obj) {
    // Call in the copy constructor
    itemsize = obj.itemsize, nbytes = obj.nbytes, stride = obj.stride,
    size_i = obj.size_i;
    dimensions = obj.dimensions;
    size = obj.size, ndim = obj.ndim;
    Sliced = obj.Sliced;

    for (int i = 0; i < dimensions; i++) {
      shape_i[i] = obj.shape_i[i];
      // cout<<"i "<<shape_i[i]<<endl;
      shape[i] = obj.shape[i];
      // cout<<"s "<<shape[i]<<endl;
      indexing[i] = obj.indexing[i];
      // cout<<"indx "<<indexing[i]<<endl;
    }

    // Copying the Data
    data = new dType[size_i];
    for (int i = 0; i < size_i; i++) {
      data[i] = obj.data[i];
    }
  }

/*
  tensor<dType, 1> operator*(tensor &obj) {
    if (shape[1] == obj.shape[0]) {
      cout<<"I am in!"<<endl;
      // Assuring if both are vectors
      if (shape[0] == 1 && obj.ndim == 1) {
        cout<<"I am in ah!"<<endl;
        int res = 0;
        for (int j = 0; j < shape[1]; j++) {
          int loc[] = {j};
          res += (*this)[loc] * obj[loc];
        }
        int re[] = {res};
        int sh[] = {1};
        tensor<dType, 1> result(re, sh);
        return result;
      } else {
        cout<<"I am in huh!"<<endl;
        int d[shape[0]];
        for (int i = 0; i < shape[0]; i++) {
          d[i] = 0;
          for (int j = 0; j < shape[1]; j++){
          int loc1[] = {i, j};
            int loc2[] = {j, i};
          d[i] += (*this)[loc1] * obj[loc2];
        } 
      }
      int shape[] = {shape[0]};
      tensor<dType, 1> result(d, shape);
      return result;
    }
  }
  else cerr << "Not Multiplicable" << endl;
}
*/
tensor operator+(const tensor &obj) {
  auto a(obj);
  if (!(this->Sliced) && !(obj.Sliced) && dimEquality(obj)) {
    auto a(obj);
    // Copying the Data
    for (int i = 0; i < size_i; i++) {
      a.data[i] = data[i] + a.data[i];
    }
    return a;
  }

  if (dimEquality(obj) && ndim <= 2) {
    if (ndim == 0) {
      for (int i = 0; i < size; i++) {
        a.data[i] = data[i] + a.data[i];
      }
    } else {
      for (int i = 0; i < size_i; i++) {
      }
    }
  }
}

tensor operator-(const tensor &obj) {
  auto a(obj);
  if (!(this->Sliced) && !(obj.Sliced) && dimEquality(obj)) {
    // Copying the Data
    for (int i = 0; i < size_i; i++) {
      a.data[i] = data[i] - a.data[i];
    }
    return a;
  }
}

bool dimEquality(tensor b) {
  if (ndim != b.ndim)
    return false;
  else
    for (int i = 0; i < ndim; i++) {
      if (shape[i] != b.shape[i])
        return false;
    }
  return true;
}
};


//Functions
   
int argmax(int* arr, int dims, int* shape, int axis=0)
{
    int max_index = 0;
    int max_value = arr[0];
    if (dims == 1)
    {
        for (int i = 1; i < shape[0]; i++)
        {
            if (arr[i] > max_value) 
            {
                max_index = i;
                max_value = arr[i];
            }
        }
    } 
    else 
    {
        int size = shape[axis];
        int strides = 1;
        for (int i = axis + 1; i < dims; i++) 
        {
            strides *= shape[i];
        }
        for (int i = 0; i < size; i++) 
        {
            int index = i * strides;
            int value = arr[index];
            for (int j = 1; j < strides; j++)
            {
              if (arr[index + j] > value)
              {
                index = index + j;
                value = arr[index];
              }
            }
            if (value > max_value)
            {
              max_index = index;
              max_value = value;
            }
        }
    }
    return max_index;
}

template<typename T, int N>
bool all(T (&arr)[N], int num_dims, int (&shape)[N], int axis) 
{
    int size = 1;
    for (int i = 0; i < num_dims; i++) {
        size *= shape[i];
    }

    if (axis < 0) {
        axis += num_dims;
    }

    int num_elements_in_axis = shape[axis];
    int num_iterations = size / num_elements_in_axis;

    for (int i = 0; i < num_iterations; i++) {
        int start_index = i * num_elements_in_axis;
        int end_index = start_index + num_elements_in_axis;

        for (int j = start_index; j < end_index; j++) {
            if (!arr[j]) {
                return false;
            }
        }
    }

    return true;
}


//numpy.argmin
void argmin(const int* arr, const int* shape, int dim, int* result, int axis =0) 
{
  int len = 1;
  for (int i = 0; i < dim; ++i) {
    len *= shape[i];
  }
  if (axis == -1) {
    int min_val = INT_MAX;
    int min_idx = 0;
    for (int i = 0; i < len; ++i) {
      if (arr[i] < min_val) {
        min_val = arr[i];
        min_idx = i;
      }
    }
    int curr_idx = dim - 1;
    int stride = 1;
    for (int i = dim - 1; i >= 0; --i) {
      int idx = min_idx / stride;
      result[curr_idx--] = idx;
      min_idx -= idx * stride;
      stride *= shape[i];
    }
    return;
  }

  int axis_len = shape[axis];
  int sub_len = len / axis_len;
  int sub_shape[dim - 1];
  int sub_stride = 1;
  for (int i = 0; i < axis; ++i) {
    sub_shape[i] = shape[i];
    sub_stride *= shape[i];
  }
  for (int i = axis + 1; i < dim; ++i) {
    sub_shape[i - 1] = shape[i];
    sub_stride *= shape[i];
  }
  int sub_arr[sub_len];
  int curr_idx = 0;
  for (int i = 0; i < axis_len; ++i) {
    memcpy(sub_arr, arr + i * sub_stride, sizeof(int) * sub_len);
    int sub_result[dim - 1];
    argmin(sub_arr, sub_shape, dim - 1, sub_result, -1);
    result[curr_idx++] = i;
    memcpy(result + curr_idx, sub_result, sizeof(int) * (dim - 1));
    curr_idx += dim - 1;
  }
}

template<typename T, size_t N>
void copy_array(const T(&input)[N], T(&output)[N]) {
  for (int i = 0; i < N; ++i) {
    output[i] = input[i];
  }
}
/*
int arr[24] = {3, 5, 2, 1, 9, 8, 6, 4, 7, 0,
                 11, 12, 10, 14, 15, 13, 19, 18, 16, 17,
                 20, 21, 22, 23};
  int shape[3] = {2, 3, 4};
  int dim = 3;
  int axis = -1;
  int result[dim];
  argmin(arr, shape, dim, axis, result);
  cout << "The minimum value is located at: ";
  for (int i = 0; i < dim; ++i) {
  cout << result[i] << " ";
  }
  cout << endl;
  return 0;

*/
template<typename T>
void fill(T *arr, const int *shape, int ndims, T scalar) {
  int size = 1;
  for (int i = 0; i < ndims; i++) {
    size *= shape[i];
  }
  for (int i = 0; i < size; i++) {
    arr[i] = scalar;
  }
}


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
