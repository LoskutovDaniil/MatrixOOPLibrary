#include "../s21_matrix_oop.h"

//-------------------------КОНСТРУКТОРЫ-------------------------//

S21Matrix::S21Matrix() : rows_(0), cols_(0), matrix_(nullptr) {}

S21Matrix::S21Matrix(int rows, int cols) {
  if (rows <= 0 || cols <= 0) {
    throw std::invalid_argument(
        "Invalid arguments: rows and cols must be positive.");
  } else {
    rows_ = rows;
    cols_ = cols;
    MemoryMatrix();
  }
}

S21Matrix::S21Matrix(const S21Matrix &other) {
  rows_ = other.rows_;
  cols_ = other.cols_;

  if (rows_ <= 0 || cols_ <= 0) {
    throw std::invalid_argument(
        "Invalid arguments: rows and cols must be positive.");
  } else {
    MemoryMatrix();

    for (int k = 0; rows_ > k; k++) {
      for (int g = 0; cols_ > g; g++) {
        matrix_[k][g] = other.matrix_[k][g];
      }
    }
  }
}

S21Matrix::S21Matrix(S21Matrix &&other) {
  rows_ = other.rows_;
  cols_ = other.cols_;
  matrix_ = other.matrix_;

  other.rows_ = 0;
  other.cols_ = 0;
  other.matrix_ = nullptr;
}

S21Matrix::~S21Matrix() {
  FreeMatrix();

  rows_ = 0, cols_ = 0;
}

//-------------------------ФУНКЦИИ-------------------------//

bool S21Matrix::EqMatrix(const S21Matrix &other) {
  if (!(SupportEqSizeMatrix(other)) || (cols_ == 0 && other.rows_ == 0))
    return false;
  else
    for (int i = 0; rows_ > i; i++) {
      for (int j = 0; cols_ > j; j++) {
        if (fabs(matrix_[i][j] - other.matrix_[i][j]) > 1e-07) return false;
      }
    }

  return true ? true : false;
}

void S21Matrix::SumMatrix(const S21Matrix &other) {
  if (!(SupportEqSizeMatrix(other))) {
    throw std::invalid_argument(
        "Invalid arguments: sum matrices are not equal.");
  }

  for (int i = 0; rows_ > i; i++) {
    for (int j = 0; cols_ > j; j++) {
      matrix_[i][j] = matrix_[i][j] + other.matrix_[i][j];
    }
  }
}

void S21Matrix::SubMatrix(const S21Matrix &other) {
  if (!(SupportEqSizeMatrix(other)))
    throw std::invalid_argument(
        "Invalid arguments: sum matrices are not equal.");
  else {
    for (int i = 0; rows_ > i; i++) {
      for (int j = 0; cols_ > j; j++) {
        matrix_[i][j] = matrix_[i][j] - other.matrix_[i][j];
      }
    }
  }
}

void S21Matrix::MulNumber(const double num) {
  for (int i = 0; rows_ > i; i++) {
    for (int j = 0; cols_ > j; j++) {
      matrix_[i][j] = matrix_[i][j] * num;
    }
  }
}

void S21Matrix::MulMatrix(const S21Matrix &other) {
  if (cols_ == other.rows_) {
    S21Matrix result(cols_, other.rows_);

    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols_; j++) {
        for (int k = 0; k < cols_; k++) {
          result.matrix_[i][j] += matrix_[i][k] * other.matrix_[k][j];
        }
      }
    }
    *this = std::move(result);
  } else {
    throw std::invalid_argument(
        "Invalid arguments: the columns of the first matrix are not equal to "
        "the rows of the second matrix");
  }
}

S21Matrix S21Matrix::Transpose() {
  S21Matrix result(cols_, rows_);
  for (int i = 0; rows_ > i; i++) {
    for (int j = 0; cols_ > j; j++) {
      result.matrix_[j][i] = matrix_[i][j];
    }
  }

  return result;
}

double S21Matrix::Determinant() {
  double result = 0.0;

  if (cols_ != rows_)
    throw std::invalid_argument("Invalid arguments: cols != rows");
  else if (rows_ == 1 && cols_ == 1)
    result = matrix_[0][0];
  else if (rows_ == 2 && cols_ == 2)
    result = matrix_[0][0] * matrix_[1][1] - matrix_[1][0] * matrix_[0][1];
  if (rows_ > 2) {
    int sign = 1;
    for (int i = 0; rows_ > i; i++) {
      S21Matrix minor = SupportSearchMinor(i, 0);
      double result_two = minor.Determinant();
      result += sign * matrix_[0][i] * result_two;
      sign *= -1;
    }
  }
  return result;
}

S21Matrix S21Matrix::CalcComplements() {
  if (rows_ != cols_)
    throw std::invalid_argument("Invalid arguments: cols != rows");

  S21Matrix result(rows_, cols_);
  {
    for (int i = 0; rows_ > i; i++) {
      for (int j = 0; cols_ > j; j++) {
        S21Matrix minor = SupportSearchMinor(i, j);
        double result_det = minor.Determinant();
        result.matrix_[j][i] = std::pow(-1, i + j) * result_det;
      }
    }
  }
  return result;
}

S21Matrix S21Matrix::InverseMatrix() {
  if (cols_ != rows_)
    throw std::invalid_argument("Invalid arguments: cols != rows");

  double result_det = Determinant();
  if (result_det == 0)
    throw std::invalid_argument("Invalid arguments: result_det = 0");

  S21Matrix calc = CalcComplements();
  S21Matrix transpose = calc.Transpose();
  transpose.MulNumber(1 / result_det);

  return transpose;
}

//-------------------------ОПЕРАТОРЫ-------------------------//

S21Matrix S21Matrix::operator+(const S21Matrix &other) {
  S21Matrix temp(*this);
  temp.SumMatrix(other);
  return temp;
}

S21Matrix S21Matrix::operator-(const S21Matrix &other) {
  S21Matrix temp(*this);
  temp.SubMatrix(other);
  return temp;
}

S21Matrix S21Matrix::operator*(const S21Matrix &other) {
  S21Matrix temp(*this);
  temp.MulMatrix(other);
  return temp;
}

S21Matrix S21Matrix::operator*(const double num) {
  S21Matrix temp(*this);
  temp.MulNumber(num);
  return temp;
}

S21Matrix &S21Matrix::operator=(const S21Matrix &other) {
  S21Matrix tmp(other);
  *this = std::move(tmp);
  return *this;
}

S21Matrix &S21Matrix::operator=(S21Matrix &&other) {
  this->~S21Matrix();
  if (&other != this) {
    rows_ = other.rows_;
    cols_ = other.cols_;
    matrix_ = other.matrix_;

    other.rows_ = 0;
    other.cols_ = 0;
    other.matrix_ = nullptr;
  }
  return *this;
}

bool S21Matrix::operator==(const S21Matrix &other) {
  return this->EqMatrix(other);
}

S21Matrix &S21Matrix::operator+=(const S21Matrix &other) {
  return *this = *this + other;
}

S21Matrix &S21Matrix::operator-=(const S21Matrix &other) {
  return *this = *this - other;
}

S21Matrix &S21Matrix::operator*=(const S21Matrix &other) {
  return *this = *this * other;
}

S21Matrix &S21Matrix::operator*=(const double num) {
  return *this = *this * num;
}

double &S21Matrix::operator()(int i, int j) {
  if (i < rows_ && 0 <= rows_ && j < cols_ && 0 <= cols_)
    return matrix_[i][j];
  else
    throw std::invalid_argument("Invalid arguments: you are out of range");
}

//-------------------------ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ-------------------------//

void S21Matrix::SetRows(int rows) {
  if (rows < 0) {
    throw std::invalid_argument("Invalid arguments: rows < 0");
  } else {
    S21Matrix tmp(rows, cols_);
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols_; j++) {
        tmp.matrix_[i][j] = matrix_[i][j];
      }
    }
    *this = tmp;
  }
}
void S21Matrix::SetCols(int cols) {
  if (cols < 0) {
    throw std::invalid_argument("Invalid arguments: cols < 0");
  } else {
    S21Matrix tmp(rows_, cols);
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols_; j++) {
        tmp.matrix_[i][j] = matrix_[i][j];
      }
    }
    *this = tmp;
  }
}

int S21Matrix::GetRows() { return rows_; }

int S21Matrix::GetCols() { return cols_; }

double **S21Matrix::GetMatrix() { return matrix_; }

void S21Matrix::FreeMatrix() {
  for (int i = 0; i < rows_; ++i) {
    delete[] matrix_[i];
  }
  delete[] matrix_;
}

void S21Matrix::MemoryMatrix() {
  matrix_ = new double *[rows_]();
  for (int i = 0; i < rows_; ++i) {
    matrix_[i] = new double[cols_]();
  }
}

S21Matrix S21Matrix::SupportSearchMinor(int pass_column, int pass_row) {
  S21Matrix minor(rows_ - 1, cols_ - 1);

  int a = 0, b = 0;

  for (int i = 0; i < rows_; i++) {
    if (i != pass_row) {
      for (int j = 0; j < cols_; j++) {
        if (j == pass_column) continue;
        minor.matrix_[a][b] = matrix_[i][j];

        if (b < minor.cols_ - 1)
          b++;
        else if (a < minor.rows_ - 1) {
          a++;
          b = 0;
        }
      }
    }
  }
  return minor;
}

bool S21Matrix::SupportEqSizeMatrix(const S21Matrix &B) {
  if ((cols_ != B.cols_ || rows_ != B.rows_)) return false;

  return true;
}