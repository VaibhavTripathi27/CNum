# CNum

In the words of wikipedia, a tensor is an algebraic object that describes a multilinear relationship between sets of algebraic objects related to a vector space. Tensors may map between different objects such as vectors, scalars, and even other tensors.

```
Basically, a tensor is just a more generalized form of a vector (or matrix), extending beyond the constraints on dimensions
```
A general example can be how the `ndarray` of *NumPy* behave.

Similar to the implementaton of `ndarray`, this project implements `tensor` thus adding support for large, multi-dimensional arrays and matrices like mathematical constructs, along with a collection of high-level mathematical functions to operate on these constructs.

`tensor` are **immutable**, **homogeneously typed** and **strided** views on memory that behave as a *n-dimensional array*, data structure.
They enable for:
- Easy operation and manipulation.
- Mathematical processing.
- Analysis and processing of large-dimensional data.
- Implementaton of mathematical models
- Pipelining processes.
- Broadcasting and Vectorization.
- All of this at the speed of C++.

## Broadcasting
The term **broadcasting** refers to behaviour of tensors with different dimensions during arithmetic operations which lead to certain constraints, the smaller array is broadcast across the larger array so that they have compatible shapes. 
![Image](/Broadcast.png "Broadcasting")


## API Reference

### Creation Routine

```cpp
  tensor<dType, dims> a(dType* data, int shape[dims]);
```

| Parameter | Type     | Description                |
| :-------- | :------- | :------------------------- |
| `dType` | `Generic Data Type` | **Required** The data type of the elements (`int`, `float`, `double`) |
| `dims` | `int` | **Required**. Number of dimensions of the tensor (should be available at compile-time) |
| `data` | `dType*` | **Required**. Pointer to a 1D array of the data to be stored in the tensor |
| `shape` | `int[dims]` | **Required**. Shape of the tensor |
| `pdim` | `int` | The dimensionality of parent tensor (for manual slicing) |
| `indexing_arr` | `int[]` | Indexing array of parent tensor |
| `srti` | `int` | Stride of the tensor |

#### Public Attributes
| Parameter | Type     | Description                |
| :-------- | :------- | :------------------------- |
| `ndim` | `int` | Dimensions of the tensor |
| `shape` | `int[ndim]` | Shape of the tensor |
| `indexing` | `int[ndim]` | Indexing coefficients |
| `size` | `int` | Total number of elements in the tensor |


#### Private Attributes


| Parameter | Type     | Description                       |
| :-------- | :------- | :-------------------------------- |
| `itemsize`      | `int` | Size of each element |
| `nbytes`      | `int` | Total size of the tensor |
| `stride`      | `int` | Stride for data access |
| `dimensions`      | `int` | Dimensionality of the parent tensor |

---

### Operators

#### Indexing of data

Takes an `int[ndim]` as input refering to location of element and returns the element at that location.

```
[0, 1, 2, 3]
[4, 5, 6, 7]
[8, 9, 10, 11]
[12, 13, 14, 15]
```
```cpp
int loc[ndim] = {2, 3};
cout << tensor_[loc];
```
```
>>11
```
#### Equality Operator `=`
Makes a deep copy of equated tensor into the invoking tensor.
```cpp
auto a = b; // a and b are tensor variables such that a is a deep copy of b
```
#### Addition or Subtraction `+` or `-`
Performs the operation and returns a tensor. If dimensions and shapes are not as per the requirements of the operation then the operand is **broadcasted** to a suitable dimensions and shape, if possible.
```cpp
auto a = b + c; // a, b and c are tensor variables such that elements a are result of element wise addition of the elements of b and c.
```
#### Asterisk `*`
Performs element-wise product of two tensors and returns a new tensor.
```cpp
auto a = b * c; // a, b and c are tensor variables such that elements a are result of element wise product of the elements of b and c.
```
---

### Functions

---

### linalg

1. dot - tensor.dot(a, b, int* result, int* shape1, int* shape2, int                     ndims1, int ndims2, int axis)
      Dot product of two arrays.
      f both a and b are 1-D arrays, it is inner product of vectors.
      If both a and b are 2-D arrays, it is matrix multiplication.
      For further range, we will stretch it's operation over n-            Dimension array.

    Input Parameters:-
    a, b - Matrix which are to be multiplied
    result - Resultant Matrix to store the resultant matrix              comprising of dot product of the given matrix
    shape1, shape2 - Matrix/Array which has the shape of the input       matrix
    ndims1,ndims2 - This attribute has the value of dimensions of        the input matrix
    Axis - This tells the axis along which multiplication has to         take place.

    Output -
    Result - Resultant matrix which has the value of the dot product     of the input matrices.

    Raises -
    Value error, if the last dimension of a is not the same size as      the second-to-last dimension of b.
   
2. vdot - vdot(vector1, vector2, length)
        Return the dot product of two vectors. vdot handles             multidimensional arrays differently than dot: it does not perform    a matrix product, but flattens input arguments to 1-D vectors        first. Consequently, it should only be used for vectors.
   Input Parameters: -
   vector1, vector2 - 2 vectors in 1D array format
   length - This parameter tells the length of the input vectors


   Output -
   result - variable containing the dot product of the vectors given.

    Further reference we will be working on the n-dimensional array      as a parameter for vectors and perform dot product of the            vectors 
3. rank - rank(a, s, dims)
    Return matrix rank of array.
    Rank of the array is the number of singular values of the array.

    Input parameters-
    a - matrix whose rank is to be found. Currently only single              matrix can be handled.
    s - shape of the matrix provided for the operation
    dims - Variable storing the dimension of the matrix passed as a             parameter

   output -
    result - Rank of the matrix

   Further reference we will be working on more than one matrix         which can be passed as parameter in the form of stack of matrices    and then we compute the rank of the matrix.
   
4. trace - trace(a, s, dims)
    Return the sum along diagonals of the array.

    Input Parameters -
    
    a - matrix whose trace is to be found.
    s - shape of the matrix provided for the operation.
    dims - Variable storing the dimension of the matrix passed as a             parameter for finding the trace.
    output
      trace - If a is 2-D, the sum along the diagonal is returned

   Further reference we will be build the trace function to be          capable of handling more than 2 dimensions, in that condition we     will be using the parameters like Axis1 and Axis2 for finding out    the 2-D sub-arrays whose traces will be returned.
   
5. solve - solve(a, b, n)

   Solve a linear matrix equation, or system of linear scalar
   equations using Gaussian Elimination method.

   Input parameters:
   
   a - Coefficient matrix.
   b - Ordinate or “dependent variable” values.

   n - It tells number of variable to solve for in the given            equations

   Output -
   Result - Solution to the system a x = b. Returned shape is           identical to b.

   Error - If a is singular or not square.

6. Inv - inv(a, n)
      Compute the inverse of a matrix using Gaussian elimination           method.

    Input parameters -
   a - matrix to be inverted
   n - sqaure matrix of dimension 'n'

   Output-
    ainv(…, M, M) ndarray or matrix.
    Inverse of the matrix a.

   Error-
      If a is not square or inversion fails.


7. Determinant - determinant(arr, n)
      Compute the determinant of an array.

   Input Parameters-
    a(…, M, M)-array_like - Input array to compute determinants for.

    n -  It tells the dimension of the array passed as a parameter.

   Output -
    Returns:
        det(…) array_like - Determinant of a.
---

