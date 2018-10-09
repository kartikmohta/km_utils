#pragma once

#include <Eigen/Core>
#include <Eigen/StdVector>

namespace km
{
namespace types
{
template <typename Scalar_t, int M, int N>
using Mat = Eigen::Matrix<Scalar_t, M, N>;

template <typename Scalar_t, int N>
using Vec = Mat<Scalar_t, N, 1>;

template <typename Scalar_t>
using MatX = Mat<Scalar_t, Eigen::Dynamic, Eigen::Dynamic>;

template <typename Scalar_t>
using VecX = Vec<Scalar_t, Eigen::Dynamic>;

template <typename Scalar_t, int N>
using SquareMat = Mat<Scalar_t, N, N>;

template <typename T>
using vector_aligned = std::vector<T, Eigen::aligned_allocator<T>>;
} // namespace types
} // namespace km
