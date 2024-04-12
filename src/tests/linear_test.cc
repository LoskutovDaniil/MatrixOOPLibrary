#include "tests.h"

TEST(TestLinear, SubtestCalcComplementsError) {
  S21Matrix mat = S21Matrix(3, 2);
  S21Matrix mat2 = S21Matrix(3, 3);
  double tmp1[9] = {1, 2, 3, 0, 4, 2, 5, 2, 1};
  double tmp2[9] = {0, 10, -20, 4, -14, 8, -8, -2, 4};
  mat.MatrixFill(tmp1);
  mat2.MatrixFill(tmp2);
  ASSERT_THROW(mat.CalcComplements(), std::invalid_argument);
}

TEST(TestLinear, SubtestDeterminant1) {
  S21Matrix mat = S21Matrix(4, 4);
  double tmp1[16] = {9, 2, 2, 9, 1, 9, 3, 4, 2, 9, 2, 2, 9, 5, 5, 9};
  mat.MatrixFill(tmp1);
  ASSERT_EQ(mat.Determinant(), 567);
}

TEST(TestLinear, SubtestDeterminant2) {
  S21Matrix mat = S21Matrix(4, 4);
  double tmp1[16] = {-9, 2, 2, 9, 1, 9, -3, 4, 2, -9, 2, 2, 9, 5, -5, 9};
  mat.MatrixFill(tmp1);
  double diff = fabs(mat.Determinant() - (double)-483);
  if (diff <= 1e-7) {
    ASSERT_EQ(1, 1);
  } else {
    ASSERT_EQ(1, 0);
  }
}

TEST(TestLinear, SubtestDeterminantError) {
  S21Matrix mat = S21Matrix(4, 3);
  double tmp1[16] = {9, 2, 2, 9, 1, 9, 3, 4, 2, 9, 2, 2, 9, 5, 5, 9};
  mat.MatrixFill(tmp1);
  ASSERT_THROW(mat.Determinant(), std::invalid_argument);
}

TEST(TestLinear, SubtestInverseMatrix) {
  S21Matrix mat = S21Matrix(3, 3);
  S21Matrix mat2 = S21Matrix(3, 3);
  double tmp1[9] = {2, 5, 7, 6, 3, 4, 5, -2, -3};
  double tmp2[9] = {1, -1, 1, -38, 41, -34, 27, -29, 24};
  mat.MatrixFill(tmp1);
  mat2.MatrixFill(tmp2);
  ASSERT_EQ(mat.InverseMatrix() == mat2, true);
}

TEST(TestLinear, SubtestInverseMatrixError1) {
  S21Matrix mat = S21Matrix(3, 2);
  S21Matrix mat2 = S21Matrix(3, 3);
  double tmp1[9] = {2, 5, 7, 6, 3, 4, 5, -2, -3};
  double tmp2[9] = {1, -1, 1, -38, 41, -34, 27, -29, 24};
  mat.MatrixFill(tmp1);
  mat2.MatrixFill(tmp2);
  ASSERT_THROW(mat.InverseMatrix(), std::invalid_argument);
}

TEST(TestLinear, SubtestInverseMatrixError2) {
  S21Matrix mat = S21Matrix(3, 3);
  double tmp1[9] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  mat.MatrixFill(tmp1);
  ASSERT_THROW(mat.InverseMatrix(), std::invalid_argument);
}
