#include "Grid.h"
#include "MathUtils.h"

template struct Grid<2>;
template struct Grid<3>;

#pragma region Grid 2D

Grid<2>::Grid(const Vector2& minCoord, const Vector2& maxCoord, const Vector2i& dim)
	: minCoord(minCoord), maxCoord(maxCoord), dimension(dim),
	vx(dim.x() + 1, dim.y()), vy(dim.x(), dim.y() + 1), tempVx(dim.x() + 1, dim.y()), tempVy(dim.x(), dim.y() + 1),
	wx(dim.x() + 1, dim.y()), wy(dim.x(), dim.y() + 1), liquidWx(dim.x() + 1, dim.y()), liquidWy(dim.x(), dim.y() + 1),
	pressure(dim), liquidPhi(dim), isValid(dim.x() + 1, dim.y() + 1), solidPhi(dim.x() + 1, dim.y() + 1) {

	Vector2 interval = maxCoord - minCoord;
	dx = Vector2(Scalar(interval.x() / dimension.x()), Scalar(interval.y() / dimension.y()));
}

Scalar Grid<2>::GetLiquidPhi(const Vector2& pos) const {
	return GetLiquidPhi(pos.x(), pos.y());
}

Scalar Grid<2>::GetLiquidPhi(const Scalar& x, const Scalar& y) const {
	if (!IsInGrid(x, y)) return 3 * dx.x();

	Scalar localX = (x - minCoord.x()) / dx.x() - Scalar(0.5);
	Scalar localY = (y - minCoord.y()) / dx.y() - Scalar(0.5);
	auto [xi, fracX] = MathUtils::GetBarycentricCoord(localX, 0, dimension.x());
	auto [yi, fracY] = MathUtils::GetBarycentricCoord(localY, 0, dimension.y());
	return MathUtils::Bilerp(liquidPhi(xi, yi), liquidPhi(xi + 1, yi), liquidPhi(xi, yi + 1), liquidPhi(xi + 1, yi + 1), fracX, fracY);
}

Scalar Grid<2>::GetSolidPhi(const Vector2& pos) const {
	return GetSolidPhi(pos.x(), pos.y());
}

Scalar Grid<2>::GetSolidPhi(const Scalar& x, const Scalar& y) const {
	if (!IsInGrid(x, y)) {
		Scalar validX = MathUtils::Clamp(x, minCoord.x() + dx.x(), maxCoord.x() - dx.x());
		Scalar validY = MathUtils::Clamp(y, minCoord.y() + dx.y(), maxCoord.y() - dx.y());
		Vector2 grad = GetSolidGradient(validX, validY);
		return GetSolidPhi(validX, validY) + grad.x() * (x - validX) + grad.y() * (y - validY);
	}

	Scalar localX = (x - minCoord.x()) / dx.x();
	Scalar localY = (y - minCoord.y()) / dx.y();
	auto [xi, fracX] = MathUtils::GetBarycentricCoord(localX, 0, dimension.x() + 1);
	auto [yi, fracY] = MathUtils::GetBarycentricCoord(localY, 0, dimension.y() + 1);
	return MathUtils::Bilerp(solidPhi(xi, yi), solidPhi(xi + 1, yi), solidPhi(xi, yi + 1), solidPhi(xi + 1, yi + 1), fracX, fracY);
}

void Grid<2>::Clear() {
	pressure.Fill(0.0);
	liquidPhi.Fill(3 * dx.x());
	vx.Fill(0.0);
	wx.Fill(0.0);
	tempVx.Fill(0.0);
	vy.Fill(0.0);
	wy.Fill(0.0);
	tempVy.Fill(0.0);
}

Vector2 Grid<2>::GetSolidGradient(const Vector2& pos) const {
	return GetSolidGradient(pos.x(), pos.y());
}

Vector2 Grid<2>::GetSolidGradient(const Scalar& x, const Scalar& y) const {
	if (!IsInGrid(x, y)) {
		Scalar validX = MathUtils::Clamp(x, minCoord.x() + dx.x(), maxCoord.x() - dx.x());
		Scalar validY = MathUtils::Clamp(y, minCoord.y() + dx.y(), maxCoord.y() - dx.y());
		return GetSolidGradient(validX, validY);
	}

	Scalar eps = 1e-4;
	Scalar left = GetSolidPhi(x - eps, y);
	Scalar right = GetSolidPhi(x + eps, y);
	Scalar bottom = GetSolidPhi(x, y - eps);
	Scalar top = GetSolidPhi(x, y + eps);
	return (Vector2(right - left, top - bottom) / (2.0 * eps)).normalized();
}

#pragma endregion

#pragma region Grid 3D

Grid<3>::Grid(const Vector3& minCoord, const Vector3& maxCoord, const Vector3i& dim)
	: minCoord(minCoord), maxCoord(maxCoord), dimension(dim),
	vx(dim.x() + 1, dim.y(), dim.z()), vy(dim.x(), dim.y() + 1, dim.z()), vz(dim.x(), dim.y(), dim.z() + 1),
	tempVx(dim.x() + 1, dim.y(), dim.z()), tempVy(dim.x(), dim.y() + 1, dim.z()), tempVz(dim.x(), dim.y(), dim.z() + 1),
	wx(dim.x() + 1, dim.y(), dim.z()), wy(dim.x(), dim.y() + 1, dim.z()), wz(dim.x(), dim.y(), dim.z() + 1),
	liquidWx(dim.x() + 1, dim.y(), dim.z()), liquidWy(dim.x(), dim.y() + 1, dim.z()), liquidWz(dim.x(), dim.y(), dim.z() + 1),
	pressure(dim), liquidPhi(dim), isValid(dim.x() + 1, dim.y() + 1, dim.z() + 1), solidPhi(dim.x() + 1, dim.y() + 1, dim.z() + 1) {

	Vector3 interval = maxCoord - minCoord;
	dx = Vector3(Scalar(interval.x() / dimension.x()), 
				 Scalar(interval.y() / dimension.y()),
				 Scalar(interval.z() / dimension.z()));
}

Scalar Grid<3>::GetLiquidPhi(const Vector3& pos) const {
	return GetLiquidPhi(pos.x(), pos.y(), pos.z());
}

Scalar Grid<3>::GetLiquidPhi(const Scalar& x, const Scalar& y, const Scalar& z) const {
	if (!IsInGrid(x, y, z)) return 3 * dx.x();

	Scalar localX = (x - minCoord.x()) / dx.x() - Scalar(0.5);
	Scalar localY = (y - minCoord.y()) / dx.y() - Scalar(0.5);
	Scalar localZ = (z - minCoord.z()) / dx.z() - Scalar(0.5);

	auto [xi, fracX] = MathUtils::GetBarycentricCoord(localX, 0, dimension.x());
	auto [yi, fracY] = MathUtils::GetBarycentricCoord(localY, 0, dimension.y());
	auto [zi, fracZ] = MathUtils::GetBarycentricCoord(localZ, 0, dimension.z());

	return MathUtils::Trilerp(
		liquidPhi(xi, yi, zi), liquidPhi(xi + 1, yi, zi),
		liquidPhi(xi, yi + 1, zi), liquidPhi(xi + 1, yi + 1, zi),
		liquidPhi(xi, yi, zi + 1), liquidPhi(xi + 1, yi, zi + 1),
		liquidPhi(xi, yi + 1, zi + 1), liquidPhi(xi + 1, yi + 1, zi + 1),
		fracX, fracY, fracZ
	);
}

Scalar Grid<3>::GetSolidPhi(const Vector3& pos) const {
	return GetSolidPhi(pos.x(), pos.y(), pos.z());
}

Scalar Grid<3>::GetSolidPhi(const Scalar& x, const Scalar& y, const Scalar& z) const {
	if (!IsInGrid(x, y, z)) {
		Scalar validX = MathUtils::Clamp(x, minCoord.x() + dx.x(), maxCoord.x() - dx.x());
		Scalar validY = MathUtils::Clamp(y, minCoord.y() + dx.y(), maxCoord.y() - dx.y());
		Scalar validZ = MathUtils::Clamp(z, minCoord.z() + dx.z(), maxCoord.z() - dx.z());
		Vector3 grad = GetSolidGradient(validX, validY, validZ);
		return GetSolidPhi(validX, validY, validZ)
			+ grad.x() * (x - validX)
			+ grad.y() * (y - validY)
			+ grad.z() * (z - validZ);
	}

	Scalar localX = (x - minCoord.x()) / dx.x();
	Scalar localY = (y - minCoord.y()) / dx.y();
	Scalar localZ = (z - minCoord.z()) / dx.z();

	auto [xi, fracX] = MathUtils::GetBarycentricCoord(localX, 0, dimension.x() + 1);
	auto [yi, fracY] = MathUtils::GetBarycentricCoord(localY, 0, dimension.y() + 1);
	auto [zi, fracZ] = MathUtils::GetBarycentricCoord(localZ, 0, dimension.z() + 1);

	return MathUtils::Trilerp(
		solidPhi(xi, yi, zi), solidPhi(xi + 1, yi, zi),
		solidPhi(xi, yi + 1, zi), solidPhi(xi + 1, yi + 1, zi),
		solidPhi(xi, yi, zi + 1), solidPhi(xi + 1, yi, zi + 1),
		solidPhi(xi, yi + 1, zi + 1), solidPhi(xi + 1, yi + 1, zi + 1),
		fracX, fracY, fracZ
	);
}

void Grid<3>::Clear() {
	pressure.Fill(0.0);
	liquidPhi.Fill(3 * dx.x());
	vx.Fill(0.0);
	wx.Fill(0.0);
	tempVx.Fill(0.0);
	vy.Fill(0.0);
	wy.Fill(0.0);
	tempVy.Fill(0.0);
	vz.Fill(0.0);
	wz.Fill(0.0);
	tempVz.Fill(0.0);
}

Vector3 Grid<3>::GetSolidGradient(const Vector3& pos) const {
	return GetSolidGradient(pos.x(), pos.y(), pos.z());
}

Vector3 Grid<3>::GetSolidGradient(const Scalar& x, const Scalar& y, const Scalar& z) const {
	if (!IsInGrid(x, y, z)) {
		Scalar validX = MathUtils::Clamp(x, minCoord.x() + dx.x(), maxCoord.x() - dx.x());
		Scalar validY = MathUtils::Clamp(y, minCoord.y() + dx.y(), maxCoord.y() - dx.y());
		Scalar validZ = MathUtils::Clamp(z, minCoord.z() + dx.z(), maxCoord.z() - dx.z());
		return GetSolidGradient(validX, validY, validZ);
	}

	Scalar eps = 1e-4;
	Scalar left = GetSolidPhi(x - eps, y, z);
	Scalar right = GetSolidPhi(x + eps, y, z);
	Scalar bottom = GetSolidPhi(x, y - eps, z);
	Scalar top = GetSolidPhi(x, y + eps, z);
	Scalar front = GetSolidPhi(x, y, z + eps);
	Scalar back = GetSolidPhi(x, y, z - eps);

	return (Vector3(right - left, top - bottom, front - back) / (2.0 * eps)).normalized();
}

#pragma endregion