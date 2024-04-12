#include "tests.h"

TEST(TestArithmetic, SubtestSum) {
  S21Matrix mat = S21Matrix(3, 3);
  S21Matrix mat2 = S21Matrix(3, 3);
  S21Matrix mat3 = S21Matrix(3, 3);
  double tmp[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
  double tmp2[16] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  double tmp3[16] = {2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17};
  mat.MatrixFill(tmp);
  mat2.MatrixFill(tmp2);
  mat3.MatrixFill(tmp3);
  mat.SumMatrix(mat2);
  ASSERT_EQ(mat == mat3, true);
}

TEST(TestArithmetic, SubtestSumError) {
  S21Matrix mat = S21Matrix(3, 3);
  S21Matrix mat2 = S21Matrix(2, 3);
  S21Matrix mat3 = S21Matrix(3, 3);
  double tmp[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
  double tmp2[16] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  double tmp3[16] = {2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17};
  mat.MatrixFill(tmp);
  mat2.MatrixFill(tmp2);
  mat3.MatrixFill(tmp3);
  ASSERT_THROW(mat.SumMatrix(mat2), std::invalid_argument);
}

TEST(TestArithmetic, SubtestSub) {
  S21Matrix mat = S21Matrix(3, 3);
  S21Matrix mat2 = S21Matrix(3, 3);
  S21Matrix mat3 = S21Matrix(3, 3);
  double tmp[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
  double tmp2[16] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  double tmp3[16] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
  mat.MatrixFill(tmp);
  mat2.MatrixFill(tmp2);
  mat3.MatrixFill(tmp3);
  mat.SubMatrix(mat2);
  ASSERT_EQ(mat == mat3, true);
}

TEST(TestArithmetic, SubtestSubError) {
  S21Matrix mat = S21Matrix(3, 3);
  S21Matrix mat2 = S21Matrix(2, 3);
  S21Matrix mat3 = S21Matrix(3, 3);
  double tmp[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
  double tmp2[16] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  double tmp3[16] = {2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17};
  mat.MatrixFill(tmp);
  mat2.MatrixFill(tmp2);
  mat3.MatrixFill(tmp3);
  ASSERT_THROW(mat.SubMatrix(mat2), std::invalid_argument);
}

TEST(TestArithmetic, SubtestMulMatrixError) {
  S21Matrix mat = S21Matrix(3, 3);
  S21Matrix mat2 = S21Matrix(2, 2);
  S21Matrix mat3 = S21Matrix(3, 3);
  S21Matrix mat4 = S21Matrix(3, 3);
  double tmp[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
  double tmp2[16] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  double tmp3[16] = {30, 36, 42, 66, 81, 96, 102, 126, 150};
  mat.MatrixFill(tmp);
  mat2.MatrixFill(tmp2);
  mat3.MatrixFill(tmp3);
  ASSERT_THROW(mat.MulMatrix(mat2), std::invalid_argument);
}
