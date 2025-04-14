#include "MathUtils.h"


Scalar MathUtils::Lerp(const Scalar& val0, const Scalar& val1, Scalar factor) {
	return val0 + MathUtils::Clamp(factor, 0.0, 1.0) * (val1 - val0);
}

Scalar MathUtils::Bilerp(const Scalar& val00, const Scalar& val10, const Scalar& val01, const Scalar& val11, Scalar factor1, Scalar factor2) {
	Scalar val0 = MathUtils::Lerp(val00, val10, factor1);
	Scalar val1 = MathUtils::Lerp(val01, val11, factor1);
	return MathUtils::Lerp(
		MathUtils::Lerp(val00, val10, factor1),
		MathUtils::Lerp(val01, val11, factor1),
		MathUtils::Clamp(factor2, 0.0, 1.0)
	);
}

Scalar MathUtils::Trilerp(const Scalar& val000, const Scalar& val100, const Scalar& val010,
	const Scalar& val110, const Scalar& val001, const Scalar& val101,
	const Scalar& val011, const Scalar& val111,
	Scalar factor1, Scalar factor2, Scalar factor3) {

	return MathUtils::Lerp(
		MathUtils::Bilerp(val000, val100, val010, val110, factor1, factor2),
		MathUtils::Bilerp(val001, val101, val011, val111, factor1, factor2),
		MathUtils::Clamp(factor3, 0.0, 1.0)
	);
}

Vector3f MathUtils::Rotate(const Vector3f& vec, const Vector3f& axis, Scalar theta) {
	Vector3f ans = vec;
	if (!MathUtils::IsSmall(theta)) {
		Scalar ct = std::cos(theta);
		Scalar st = std::sin(theta);
		ans = ct * vec + st * axis.cross(vec) + axis.dot(vec) * (1.0 - ct) * axis;
	}
	return ans;
}

Scalar MathUtils::ComputeAngle(const Vector3f& vec0, const Vector3f& vec1, const Vector3f& axis) {
	Vector3f vec = vec0.cross(vec1);
	Scalar theta = std::atan2(vec.norm(), vec0.dot(vec1));
	if (axis.dot(vec) < 0.0) theta = -theta;
	return theta;
}

Matrix3f MathUtils::CrossMatrix(const Vector3f& vec) {
	Matrix3f mat;
	mat <<	0.0,		-vec.z(),	vec.y(),
			vec.z(),	0.0,		-vec.x(),
			-vec.y(),	vec.x(),	0.0;
	return mat;
}

Scalar MathUtils::FractionInRegion(const Scalar & left, const Scalar& right) {
	Scalar ans = 0.0;
	if (left < 0.0 && right < 0.0) ans = 1.0;
	else if (left < 0.0 && 0.0 < right) ans = left / (left - right);
	else if (left > 0.0 && 0.0 > right) ans = right / (right - left);
	return ans;
}

Scalar MathUtils::FractionInRegion(const Scalar& val00, const Scalar& val10, const Scalar& val01, const Scalar& val11) {
	if (val00 < 0.0 && val10 < 0.0 && val01 < 0.0 && val11 < 0.0) return 1.0;
	if (val00 >= 0.0 && val10 >= 0.0 && val01 >= 0.0 && val11 >= 0.0) return 0.0;
	
	std::vector<Vector2> points;
	if (val00 < 0.0) points.push_back(Vector2(0.0, 0.0));
	if (val01 < 0.0) points.push_back(Vector2(0.0, 1.0));
	if (val10 < 0.0) points.push_back(Vector2(1.0, 0.0));
	if (val11 < 0.0) points.push_back(Vector2(1.0, 1.0));

	if ((val00 < 0.0 && val10 >= 0.0) || (val00 >= 0.0 && val10 < 0.0)) {
		points.push_back(Vector2(val00 / (val00 - val10), 0.0));
	}
	if ((val00 < 0.0 && val01 >= 0.0) || (val00 >= 0.0 && val01 < 0.0)) {
		points.push_back(Vector2(0.0, val00 / (val00 - val01)));
	}
	if ((val11 < 0.0 && val01 >= 0.0) || (val11 >= 0.0 && val01 < 0.0)) {
		points.push_back(Vector2(val01 / (val01 - val11), 1.0));
	}
	if ((val11 < 0.0 && val10 >= 0.0) || (val11 >= 0.0 && val10 < 0.0)) {
		points.push_back(Vector2(1.0, val10 / (val10 - val11)));
	}

#ifdef _DEBUG
	ASSERT_MSG(!points.empty(), "points cannot be empty here!");
#endif

	Vector2 center = Vector2::Zero();
	for (const Vector2& pos : points) center += pos;
	center /= Scalar(points.size());

	std::sort(points.begin(), points.end(), [center](const Vector2& p1, const Vector2& p2) {
		return std::atan2(p1.y() - center.y(), p1.x() - center.x()) < std::atan2(p2.y() - center.y(), p2.x() - center.x());
	});

	return ComputePolygonArea(points);
}

std::pair<Integer, Scalar> MathUtils::GetBarycentricCoord(const Scalar& x, const Integer& minIndex, const Integer& maxIndex) {
	Integer index = static_cast<Integer>(std::floor(x));
	Scalar frac = x - index;

	if (index < minIndex) {
		index = minIndex;
		frac = 0.0;
	}
	else if (index > maxIndex - 2) {
		index = maxIndex - 2;
		frac = 1.0;
	}
	return { index, frac };
}

Scalar MathUtils::ComputePolygonArea(const std::vector<Vector2>& points) {
	// shoelace formula
	Scalar area = 0.0;
	Integer size = points.size();
	for (Integer index = 0; index < size; ++index) {
		const Vector2& p1 = points[index];
		const Vector2& p2 = points[(index + 1) % size];
		area += p1.x() * p2.y() - p2.x() * p1.y();
	}
	return std::abs(area) * 0.5;
}
