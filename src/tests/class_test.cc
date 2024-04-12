#include "tests.h"

TEST(TestClass, SubtestInitArgs) {
  S21Matrix mat = S21Matrix(4, 4);
  double tmp[16] = {9, 2, 2, 9, 1, 9, 3, 4, 2, 9, 2, 2, 9, 5, 5, 9};
  ASSERT_EQ(mat.GetRows(), 4);
  ASSERT_EQ(mat.GetRows(), 4);
  mat.MatrixFill(tmp);
  for (int i = 0; i < mat.GetRows(); i++) {
    for (int j = 0; j < mat.GetRows(); j++) {
      ASSERT_EQ(mat(i, j), tmp[i * mat.GetRows() + j]);
    }
  }
}

TEST(TestClass, SubtestInitOutOfRange) {
  ASSERT_THROW(S21Matrix(-1, 4), std::invalid_argument);
}

TEST(TestClass, SubtestCopy) {
  S21Matrix mat = S21Matrix(4, 4);
  double tmp[16] = {9, 2, 2, 9, 1, 9, 3, 4, 2, 9, 2, 2, 9, 5, 5, 9};
  mat.MatrixFill(tmp);
  S21Matrix mat2 = S21Matrix(mat);
  ASSERT_EQ(mat.EqMatrix(mat2), true);
}

TEST(TestClass, SubtestMove) {
  S21Matrix mat = S21Matrix(4, 4);
  double tmp[16] = {9, 2, 2, 9, 1, 9, 3, 4, 2, 9, 2, 2, 9, 5, 5, 9};
  mat.MatrixFill(tmp);
  S21Matrix mat2 = S21Matrix(std::move(mat));
  ASSERT_EQ(mat2.GetRows(), 4);
  ASSERT_EQ(mat2.GetRows(), 4);
  for (int i = 0; i < mat.GetRows(); i++) {
    for (int j = 0; j < mat.GetRows(); j++) {
      ASSERT_EQ(mat2(i, j), tmp[i * mat.GetRows() + j]);
    }
  }
  ASSERT_EQ(mat.GetRows(), 0);
  ASSERT_EQ(mat.GetRows(), 0);
}

TEST(TestClass, SubtestgetRow) {
  S21Matrix mat = S21Matrix(4, 3);
  ASSERT_EQ(mat.GetRows(), 4);
}

TEST(TestClass, SubtestgetCol) {
  S21Matrix mat = S21Matrix(3, 4);
  ASSERT_EQ(mat.GetCols(), 4);
}

TEST(TestClass, SubtestsetRow) {
  S21Matrix mat = S21Matrix(3, 4);
  S21Matrix mat2 = S21Matrix(4, 4);
  double tmp[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
  double tmp2[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
  mat2.MatrixFill(tmp2);
  mat.SetRows(4);
  mat.MatrixFill(tmp);
  ASSERT_EQ(mat.EqMatrix(mat2), true);
}

TEST(TestClass, SubtestsetRowError) {
  S21Matrix mat = S21Matrix(3, 4);
  ASSERT_THROW(mat.SetRows(-1), std::invalid_argument);
}

TEST(TestClass, SubtestsetCol) {
  S21Matrix mat = S21Matrix(4, 3);
  S21Matrix mat2 = S21Matrix(4, 4);
  double tmp[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
  double tmp2[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
  mat.SetCols(4);
  mat.MatrixFill(tmp);
  mat2.MatrixFill(tmp2);

  ASSERT_EQ(mat.EqMatrix(mat2), true);
}

TEST(TestClass, SubtestsetColError) {
  S21Matrix mat = S21Matrix(3, 4);
  ASSERT_THROW(mat.SetCols(-1), std::invalid_argument);
}

TEST(TestClass, SubtestGetValue) {
  S21Matrix mat = S21Matrix(3, 4);
  double tmp[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
  mat.MatrixFill(tmp);
  ASSERT_EQ(mat(0, 0), 1);
}

TEST(TestClass, SubtestGetValueError) {
  S21Matrix mat = S21Matrix(3, 4);
  double tmp[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
  mat.MatrixFill(tmp);
  ASSERT_THROW(mat(10, 3), std::invalid_argument);
}

TEST(TestClass, SubtestSetValue) {
  S21Matrix mat = S21Matrix(3, 4);
  double tmp[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
  mat.MatrixFill(tmp);
  mat(0, 0) = 0;
  ASSERT_EQ(mat(0, 0), 0);
}

TEST(TestClass, SubtestSetValueError) {
  S21Matrix mat = S21Matrix(3, 4);
  double tmp[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
  mat.MatrixFill(tmp);
  ASSERT_THROW(mat(10, 3) = 0, std::invalid_argument);
}

TEST(TestClass, SubtestSetCopyMatrix) {
  S21Matrix mat = S21Matrix(3, 4);
  double tmp[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
  mat.MatrixFill(tmp);
  S21Matrix mat2 = mat;
  ASSERT_EQ(mat.EqMatrix(mat2), true);
}

TEST(TestClass, SubtestSetMatrix) {
  S21Matrix mat = S21Matrix(4, 4);
  double tmp[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
  mat.MatrixFill(tmp);
  S21Matrix mat2 = mat;
  ASSERT_EQ(mat.EqMatrix(mat2), true);
}

TEST(TestClass, SubtestMoveMatrix) {
  S21Matrix mat = S21Matrix(3, 4);
  double tmp[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
  mat.MatrixFill(tmp);
  S21Matrix mat2 = std::move(mat);
  ASSERT_EQ(mat2.GetRows(), 3);
  ASSERT_EQ(mat2.GetCols(), 4);
  for (int i = 0; i < mat.GetRows(); i++) {
    for (int j = 0; j < mat.GetCols(); j++) {
      ASSERT_EQ(mat2(i, j), tmp[i * mat.GetRows() + j]);
    }
  }
  ASSERT_EQ(mat.GetRows(), 0);
  ASSERT_EQ(mat.GetRows(), 0);
}

TEST(TestClass, SubtestSum) {
  S21Matrix mat = S21Matrix(3, 3);
  S21Matrix mat2 = S21Matrix(3, 3);
  S21Matrix mat3 = S21Matrix(3, 3);
  S21Matrix mat4 = S21Matrix(3, 3);
  double tmp[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
  double tmp2[16] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  double tmp3[16] = {2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17};
  mat.MatrixFill(tmp);
  mat2.MatrixFill(tmp2);
  mat3.MatrixFill(tmp3);
  mat4 = mat + mat2;
  ASSERT_EQ(mat4.EqMatrix(mat3), true);
}

TEST(TestClass, SubtestSub) {
  S21Matrix mat = S21Matrix(4, 4);
  S21Matrix mat2 = S21Matrix(4, 4);
  S21Matrix mat3 = S21Matrix(4, 4);
  S21Matrix mat4 = S21Matrix(4, 4);
  double tmp[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
  double tmp2[16] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  double tmp3[16] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
  mat.MatrixFill(tmp);
  mat2.MatrixFill(tmp2);
  mat3.MatrixFill(tmp3);
  mat4 = mat - mat2;
  ASSERT_EQ(mat4.EqMatrix(mat3), true);
}

TEST(TestClass, SubtestMulNum) {
  S21Matrix mat = S21Matrix(3, 3);
  S21Matrix mat2 = S21Matrix(3, 3);
  S21Matrix mat3 = S21Matrix(3, 3);
  S21Matrix mat4 = S21Matrix(3, 3);
  double tmp[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
  double tmp2[16] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  double tmp3[16] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  mat.MatrixFill(tmp);
  mat2.MatrixFill(tmp2);
  mat3.MatrixFill(tmp3);
  mat4 = mat * 0;
  ASSERT_EQ(mat4.EqMatrix(mat3), true);
}

TEST(TestClass, SubtestEqTrue) {
  S21Matrix mat = S21Matrix(3, 3);
  S21Matrix mat2 = S21Matrix(3, 3);
  double tmp[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
  mat.MatrixFill(tmp);
  mat2.MatrixFill(tmp);
  ASSERT_EQ(mat == mat2, true);
}

TEST(TestClass, SubtestEqFalse) {
  S21Matrix mat = S21Matrix(3, 3);
  S21Matrix mat2 = S21Matrix(3, 3);
  double tmp[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
  mat.MatrixFill(tmp);
  mat2.MatrixFill(tmp);
  mat2(0, 0) = 0;
  ASSERT_EQ(mat == mat2, false);
}

TEST(TestClass, SubtestSumEq) {
  S21Matrix mat = S21Matrix(3, 3);
  S21Matrix mat2 = S21Matrix(3, 3);
  S21Matrix mat3 = S21Matrix(3, 3);
  double tmp[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
  double tmp2[16] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  double tmp3[16] = {2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17};
  mat.MatrixFill(tmp);
  mat2.MatrixFill(tmp2);
  mat3.MatrixFill(tmp3);
  mat += mat2;
  ASSERT_EQ(mat.EqMatrix(mat3), true);
}

TEST(TestClass, SubtestSubEq) {
  S21Matrix mat = S21Matrix(3, 3);
  S21Matrix mat2 = S21Matrix(3, 3);
  S21Matrix mat3 = S21Matrix(3, 3);
  double tmp[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
  double tmp2[16] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  double tmp3[16] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
  mat.MatrixFill(tmp);
  mat2.MatrixFill(tmp2);
  mat3.MatrixFill(tmp3);
  mat -= mat2;
  ASSERT_EQ(mat.EqMatrix(mat3), true);
}

TEST(TestClass, SubtestMulNumEq) {
  S21Matrix mat = S21Matrix(3, 3);
  S21Matrix mat2 = S21Matrix(3, 3);
  S21Matrix mat3 = S21Matrix(3, 3);
  double tmp[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
  double tmp2[16] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  double tmp3[16] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  mat.MatrixFill(tmp);
  mat2.MatrixFill(tmp2);
  mat3.MatrixFill(tmp3);
  mat *= 0;
  ASSERT_EQ(mat.EqMatrix(mat3), true);
}

TEST(TestClass, SubtestMulMatEq) {
  S21Matrix mat = S21Matrix(3, 3);
  S21Matrix mat2 = S21Matrix(3, 3);
  S21Matrix mat3 = S21Matrix(3, 3);
  double tmp[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
  double tmp2[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
  double tmp3[16] = {30, 36, 42, 66, 81, 96, 102, 126, 150};
  mat.MatrixFill(tmp);
  mat2.MatrixFill(tmp2);
  mat3.MatrixFill(tmp3);
  mat *= mat;
  ASSERT_EQ(mat.EqMatrix(mat3), true);
}
