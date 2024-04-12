#include "../s21_matrix_oop.h"

void S21Matrix::MatrixFill(double* array) {
  int k = 0;
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = array[k];
      k++;
    }
  }
}

void SupportCoutMatrix(S21Matrix& matrix) {
  int rows = matrix.GetRows();
  int cols = matrix.GetCols();
  double** matrix_ = matrix.GetMatrix();

  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      std::cout << matrix_[i][j] << " ";
    }
    std::cout << "\n";
  }
  std::cout << "\n";
}