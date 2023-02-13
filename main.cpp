#include "numC.h"

int main() {
  int data[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17};
  int shape[] = {2, 3, 3};
  tensor<int, 3> ten(data, shape);
  int data1[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18};
  tensor<int, 3> nten(data1, shape);
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        int loc[] = {i, j, k};
        cout << ten[loc] << " ";
      }
      cout << endl;
    }
    cout << endl;
    cout << endl;
  }
  
  int slice[3][2] = {{0, 0}, {1, 2}, {1, 2}};
  ten.slice(slice);
  for (int j = 0; j < 2; j++) {
    for (int k = 0; k < 2; k++) {
      int loc[] = {j, k};
      cout << ten[loc] << " ";
    }
    cout << endl;
  }

  ten.resetSlice();
  cout << ten.Sliced << " " << nten.Sliced << endl;

  auto tco = nten + ten;
  cout << tco.shape[0] << " " << tco.shape[1] << " " << tco.ndim << endl
       << endl;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        int loc[] = {i, j, k};
        cout << tco[loc] << " ";
      }
      cout << endl;
    }
    cout << endl;
    cout << endl;
  }
}
