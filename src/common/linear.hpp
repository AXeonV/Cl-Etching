#pragma once
#include <algorithm>
#include <cmath>

namespace etch {
// 对称特征分解构造 Moore-Penrose 最小二乘解，同时返回局部拟合矩阵的秩。
template <int N>
int least_squares(double (&matrix)[N][N], const double (&rhs)[N], double (&result)[N]) {
  double eigenvectors[N][N] = {};
  for (int i = 0; i < N; ++i)
    eigenvectors[i][i] = 1;
  // 每轮消去最大的非对角项，并累积对应的特征向量旋转。
  for (int iteration = 0; iteration < 40 * N * N / 9; ++iteration) {
    int p = 0, q = 1;
    for (int i = 0; i < N; ++i)
      for (int j = i + 1; j < N; ++j)
        if (std::abs(matrix[i][j]) > std::abs(matrix[p][q])) {
          p = i;
          q = j;
        }
    double scale = 0;
    for (int i = 0; i < N; ++i)
      scale += std::abs(matrix[i][i]);
    if (std::abs(matrix[p][q]) < 1e-14 * scale)
      break;
    double angle = 0.5 * std::atan2(2 * matrix[p][q], matrix[q][q] - matrix[p][p]);
    double c = std::cos(angle), s = std::sin(angle);
    double pp = matrix[p][p], qq = matrix[q][q], pq = matrix[p][q];
    for (int k = 0; k < N; ++k)
      if (k != p && k != q) {
        double kp = matrix[k][p], kq = matrix[k][q];
        matrix[k][p] = matrix[p][k] = c * kp - s * kq;
        matrix[k][q] = matrix[q][k] = s * kp + c * kq;
      }
    matrix[p][p] = c * c * pp - 2 * c * s * pq + s * s * qq;
    matrix[q][q] = s * s * pp + 2 * c * s * pq + c * c * qq;
    matrix[p][q] = matrix[q][p] = 0;
    for (int k = 0; k < N; ++k) {
      double vp = eigenvectors[k][p], vq = eigenvectors[k][q];
      eigenvectors[k][p] = c * vp - s * vq;
      eigenvectors[k][q] = s * vp + c * vq;
    }
  }
  double maximum = matrix[0][0];
  for (int i = 1; i < N; ++i)
    maximum = std::max(maximum, matrix[i][i]);

  // 只对有效特征值构造伪逆，使秩不足的拟合仍返回最小范数解。
  int rank = 0;
  for (int j = 0; j < N; ++j)
    if (matrix[j][j] > maximum * 1e-12) {
      ++rank;
      double projected = 0;
      for (int k = 0; k < N; ++k)
        projected += eigenvectors[k][j] * rhs[k];
      for (int k = 0; k < N; ++k)
        result[k] += eigenvectors[k][j] * projected / matrix[j][j];
    }
  return rank;
}

}
