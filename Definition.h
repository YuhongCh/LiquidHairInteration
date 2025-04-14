#pragma once

#define DIMENSION 3

#ifdef _DEBUG
#include <cassert>

#define ASSERT_MSG(condition, message)                   \
    do {                                                 \
        if (!(condition)) {                              \
            std::cout << message <<std::endl;            \
            assert(condition);                           \
        }                                                \
    } while (false)

#endif

#pragma warning( disable : 4244 )
#pragma warning( disable : 4305 ) 

#pragma region STL and Basics
#include <windows.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <cstddef>
#include <random>

using Scalar = float;
using Integer = int;
using UInteger = unsigned int;
#pragma endregion

#pragma region Eigen
#include <Eigen/Core>
#include <Eigen/Geometry>

using Vector2 = Eigen::Matrix<Scalar, 2, 1>;
using Vector3 = Eigen::Matrix<Scalar, 3, 1>;
using Vector4 = Eigen::Matrix<Scalar, 4, 1>;

using Vector2f = Vector2;
using Vector3f = Vector3;
using Vector4f = Vector4;

using Vector2i = Eigen::Matrix<Integer, 2, 1>;
using Vector3i = Eigen::Matrix<Integer, 3, 1>;
using Vector4i = Eigen::Matrix<Integer, 4, 1>;

// Templated aliases for column vectors:
template <int N>
using VectorXf = Eigen::Matrix<Scalar, N, 1>;
template <int N>
using VectorXi = Eigen::Matrix<Integer, N, 1>;
template <int N>
using VectorX = Eigen::Matrix<Scalar, N, 1>;

using Vector2s = std::vector<Vector2>;
using Vector3s = std::vector<Vector3>;
using Vector4s = std::vector<Vector4>;

template <int N>
using VectorXs = std::vector<VectorX<N>>;

using Matrix2 = Eigen::Matrix<Scalar, 2, 2>;
using Matrix3 = Eigen::Matrix<Scalar, 3, 3>;
using Matrix4 = Eigen::Matrix<Scalar, 4, 4>;

using Matrix2f = Matrix2;
using Matrix3f = Matrix3;
using Matrix4f = Matrix4;

using Matrix2i = Eigen::Matrix<Integer, 2, 2>;
using Matrix3i = Eigen::Matrix<Integer, 3, 3>;
using Matrix4i = Eigen::Matrix<Integer, 4, 4>;

template <int N, int M>
using MatrixXf = Eigen::Matrix<Scalar, N, M>;
template <int N, int M>
using MatrixXi = Eigen::Matrix<Integer, N, M>;
template <int N, int M>
using MatrixX = Eigen::Matrix<Scalar, N, M>;

using Matrix2s = std::vector<Matrix2>;
using Matrix3s = std::vector<Matrix3>;
using Matrix4s = std::vector<Matrix4>;

template <int N, int M>
using MatrixXs = std::vector<MatrixX<N, M>>;
#pragma endregion

#pragma region Constant
constexpr Scalar EPSILON = 1e-8;
constexpr Scalar PI = 3.1415926535897932;
#pragma endregion

#pragma region HPC
#include <omp.h>

#pragma endregion

#pragma region Common Utility
#define BUFFER_OFFSET(offset) ((void *) (offset))

inline Integer RandomRange(Integer minValue, Integer maxValue) { return std::rand() % (maxValue - minValue) + minValue; }
inline Scalar RandomRange(Scalar minValue, Scalar maxValue) { return minValue + static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX) * (maxValue - minValue); }

template <typename Derived>
std::string Eigen2String(const Eigen::MatrixBase<Derived>& vec) {
	std::stringstream ss;
	ss << vec;
	return ss.str();
}

#pragma endregion

#pragma region Setting


#pragma endregion