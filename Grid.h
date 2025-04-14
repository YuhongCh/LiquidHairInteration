#pragma once

#include "Definition.h"


template <int DIM>
struct Grid;

#pragma region Grid 2D
template <>
struct Grid<2> {
public:

	Grid(const Vector2& minCoord, const Vector2& maxCoord, const Vector2i& dimension);

	inline bool IsInGrid(const Vector2& pos) const {
		return minCoord.x() <= pos.x() && minCoord.y() <= pos.y() && pos.x() <= maxCoord.x() && pos.y() <= maxCoord.y();
	}

	inline bool IsInGrid(const Scalar& x, const Scalar& y) const {
		return minCoord.x() <= x && minCoord.y() <= y && x <= maxCoord.x() && y <= maxCoord.y();
	}

	inline bool IsValidIndex(const Vector2i& index) const {
		return 0 <= index.x() && 0 <= index.y() && index.x() < dimension.x() && index.y() < dimension.y();
	}

	inline bool IsValidIndex(Integer xi, Integer yi) const {
		return 0 <= xi && 0 <= yi && xi < dimension.x() && yi < dimension.y();
	}

	Scalar GetLiquidPhi(const Vector2& pos) const;

	Scalar GetLiquidPhi(const Scalar& x, const Scalar& y) const;

	Scalar GetSolidPhi(const Vector2& pos) const;

	Scalar GetSolidPhi(const Scalar& x, const Scalar& y) const;

	inline Vector2 Index2Point(const Vector2i& index) const { return (index.cast<Scalar>() + Vector2(0.5, 0.5)).cwiseProduct(dx) + minCoord; }

	inline Vector2 Index2Point(const Integer& xi, const Integer& yi) const { return Index2Point(Vector2i(xi, yi)); }

	inline Vector2i Point2Index(const Vector2& pos) const { return Vector2i(
			static_cast<Integer>((pos.x() - minCoord.x()) / dx.x()), 
			static_cast<Integer>((pos.y() - minCoord.y()) / dx.y())
		); 
	}

	inline Vector2i Point2Index(const Scalar& x, const Scalar& y) const { return Point2Index(Vector2(x, y)); }

	void Clear();

	Vector2 GetSolidGradient(const Vector2& pos) const;

	Vector2 GetSolidGradient(const Scalar& x, const Scalar& y) const;

protected:

public:
	Vector2 minCoord;
	Vector2 maxCoord;
	Vector2i dimension;

	Vector2 dx;
	
	std::vector<std::vector<Scalar>> vx;
	std::vector<std::vector<Scalar>> vy;

	// for extrapolate
	std::vector<std::vector<Scalar>> tempVx;
	std::vector<std::vector<Scalar>> tempVy;
	std::vector<std::vector<bool>> isValid;

	std::vector<std::vector<Scalar>> wx;
	std::vector<std::vector<Scalar>> wy;
	std::vector<std::vector<Scalar>> liquidWx;
	std::vector<std::vector<Scalar>> liquidWy;

	std::vector<std::vector<Scalar>> pressure;
	std::vector<std::vector<Scalar>> liquidPhi;
	std::vector<std::vector<Scalar>> solidPhi;
	
};

#pragma endregion

#pragma region Grid 3D
template <>
struct Grid<3> {
public:

	Grid(const Vector3& minCoord, const Vector3& maxCoord, const Vector3i& dimension);

	inline bool IsInGrid(const Vector3& pos) const {
		return minCoord.x() <= pos.x() && minCoord.y() <= pos.y() && minCoord.z() <= pos.z() &&
			   pos.x() <= maxCoord.x() && pos.y() <= maxCoord.y() && pos.z() <= maxCoord.z();
	}

	inline bool IsInGrid(const Scalar& x, const Scalar& y, const Scalar& z) const {
		return minCoord.x() <= x && minCoord.y() <= y && minCoord.z() <= z &&
			   x <= maxCoord.x() && y <= maxCoord.y() && z <= maxCoord.z();
	}

	inline bool IsValidIndex(const Vector3i& index) const {
		return 0 <= index.x() && 0 <= index.y() && 0 <= index.z() && 
			   index.x() < dimension.x() && index.y() < dimension.y() && index.z() < dimension.z();
	}

	inline bool IsValidIndex(Integer xi, Integer yi, Integer zi) const {
		return 0 <= xi && 0 <= yi && 0 <= zi && xi < dimension.x() && yi < dimension.y() && zi < dimension.z();
	}

	Scalar GetLiquidPhi(const Vector3& pos) const;

	Scalar GetLiquidPhi(const Scalar& x, const Scalar& y, const Scalar& z) const;

	Scalar GetSolidPhi(const Vector3& pos) const;

	Scalar GetSolidPhi(const Scalar& x, const Scalar& y, const Scalar& z) const;

	inline Vector3 Index2Point(const Vector3i& index) const { return (index.cast<Scalar>() + Vector3(0.5, 0.5, 0.5)).cwiseProduct(dx) + minCoord; }

	inline Vector3 Index2Point(const Integer& xi, const Integer& yi, const Integer& zi) const { return Index2Point(Vector3i(xi, yi, zi)); }

	inline Vector3i Point2Index(const Vector3& pos) const {
		return Vector3i(
			static_cast<Integer>((pos.x() - minCoord.x()) / dx.x()),
			static_cast<Integer>((pos.y() - minCoord.y()) / dx.y()),
			static_cast<Integer>((pos.z() - minCoord.z()) / dx.z())
		);
	}

	inline Vector3i Point2Index(const Scalar& x, const Scalar& y, const Scalar& z) const { return Point2Index(Vector3(x, y, z)); }

	void Clear();

	Vector3 GetSolidGradient(const Vector3& pos) const;

	Vector3 GetSolidGradient(const Scalar& x, const Scalar& y, const Scalar& z) const;

protected:

public:
	Vector3 minCoord;
	Vector3 maxCoord;
	Vector3i dimension;

	Vector3 dx;

	std::vector<std::vector<std::vector<Scalar>>> vx;
	std::vector<std::vector<std::vector<Scalar>>> vy;
	std::vector<std::vector<std::vector<Scalar>>> vz;

	// for extrapolate
	std::vector<std::vector<std::vector<Scalar>>> tempVx;
	std::vector<std::vector<std::vector<Scalar>>> tempVy;
	std::vector<std::vector<std::vector<Scalar>>> tempVz;
	std::vector<std::vector<std::vector<bool>>> isValid;

	std::vector<std::vector<std::vector<Scalar>>> wx;
	std::vector<std::vector<std::vector<Scalar>>> wy;
	std::vector<std::vector<std::vector<Scalar>>> wz;
	std::vector<std::vector<std::vector<Scalar>>> liquidWx;
	std::vector<std::vector<std::vector<Scalar>>> liquidWy;
	std::vector<std::vector<std::vector<Scalar>>> liquidWz;

	std::vector<std::vector<std::vector<Scalar>>> pressure;
	std::vector<std::vector<std::vector<Scalar>>> liquidPhi;
	std::vector<std::vector<std::vector<Scalar>>> solidPhi;

};

#pragma endregion