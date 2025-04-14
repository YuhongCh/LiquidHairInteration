#pragma once
#include "Definition.h"


class MathUtils {
public:

	inline static bool IsSmall(const Scalar& val) { return std::abs(val) < EPSILON; }

	inline static Scalar Clamp(const Scalar& val, const Scalar& minVal, const Scalar& maxVal) { return min(maxVal, max(minVal, val)); }

	static Scalar Lerp(const Scalar& val0, const Scalar& val1, Scalar factor);

	static Scalar Bilerp(const Scalar& val00, const Scalar& val10, const Scalar& val01, const Scalar& val11, Scalar factor1, Scalar factor2);

	static Scalar Trilerp(const Scalar& val000, const Scalar& val100, const Scalar& val010,
					const Scalar& val110, const Scalar& val001, const Scalar& val101,
					const Scalar& val011, const Scalar& val111,
					Scalar factor1, Scalar factor2, Scalar factor3);


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

};
