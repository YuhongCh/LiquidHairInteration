#pragma once
#include "Definition.h"


class MathUtils {
public:

	inline static bool IsSmall(const Scalar& val) { return std::abs(val) < EPSILON; }

	inline static Scalar Clamp(const Scalar& val, const Scalar& minVal, const Scalar& maxVal) { return min(maxVal, max(minVal, val)); }

	template <typename T>
	static T Lerp(const T& val0, const T& val1, Scalar factor) {
		return val0 + Clamp(factor, 0.0, 1.0) * (val1 - val0);
	}

	template <typename T>
	static T Bilerp(const T& val00, const T& val10, const T& val01, const T& val11, Scalar factor1, Scalar factor2) {
		T val0 = Lerp(val00, val10, factor1);
		T val1 = Lerp(val01, val11, factor1);
		return Lerp(val0, val1, Clamp(factor2, 0.0, 1.0));
	}

	template <typename T>
	static T Trilerp(const T& val000, const T& val100, const T& val010,
					 const T& val110, const T& val001, const T& val101,
					 const T& val011, const T& val111,
					 Scalar factor1, Scalar factor2, Scalar factor3) {
		return Lerp(
			Bilerp(val000, val100, val010, val110, factor1, factor2),
			Bilerp(val001, val101, val011, val111, factor1, factor2),
			Clamp(factor3, 0.0, 1.0)
		);
	}


	template <int N>
	inline static VectorXf<N> OrthoNormalized(const VectorXf<N>& vec0, const VectorXf<N>& vec1) { return (vec0 - vec0.dot(vec1) * vec1).normalized(); }

	static Vector3f Rotate(const Vector3f& vec, const Vector3f& axis, Scalar theta);

	static Scalar ComputeAngle(const Vector3f& vec0, const Vector3f& vec1, const Vector3f& axis);

	static Matrix3f CrossMatrix(const Vector3f& vec);

	static Scalar FractionInRegion(const Scalar& left, const Scalar& right);

	static Scalar FractionInRegion(const Scalar& val00, const Scalar& val10, const Scalar& val01, const Scalar& val11);

	static std::pair<Integer, Scalar> GetBarycentricCoord(const Scalar& x, const Integer& minIndex, const Integer& maxIndex);
	
	template <typename T>
	static inline T reduce(const std::vector<T>& vec0, const std::vector<T>& vec1) {
		assert(vec0.size() == vec1.size());
		T sum = 0;

		#pragma omp parallel for reduction(+:sum)
		for (size_t i = 0; i < vec0.size(); i++) {
			sum += vec0[i] * vec1[i];
		}
		return sum;
	}

	// MUST enfore points's order in all CCW or CW
	static Scalar ComputePolygonArea(const std::vector<Vector2>& points);

	template <int DIM>
	static VectorX<DIM> ComputeNormal(const VectorX<DIM>& vec);

	/// <summary>
	/// Parallel transport prevFrame to currFrame based on prevTangent and currTangent
	/// IMPORTANT: assume all vectors are normalized
	/// </summary>
	template <int DIM>
	static VectorX<DIM> ParallelTransport(const VectorX<DIM>& prevFrame, const VectorX<DIM>& prevTangent, const VectorX<DIM>& currTangent);

};
