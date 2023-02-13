#include "main.cpp"
#include <vector>
#include <algorithm>
#include <cmath>
#include <climits>
#include <cstring>

//Member functions for ndarray


   
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

//numpy.all

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
    argmin(sub_arr, sub_shape, dim - 1, -1, sub_result);
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

//numpy.ndarray.fill:

void fill(T *arr, const int *shape, int ndims, T scalar) {
  int size = 1;
  for (int i = 0; i < ndims; i++) {
    size *= shape[i];
  }
  for (int i = 0; i < size; i++) {
    arr[i] = scalar;
  }
}


