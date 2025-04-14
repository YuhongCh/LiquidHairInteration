#pragma once

#include "Grid.h"


#pragma region Jacobi Solver
template <int DIM>
class JacobiSolver;

#pragma region Jacobi Solver 2D
template <>
class JacobiSolver<2> {
public:
	JacobiSolver(Grid<2>& grid, Integer maxIteration = 500, const Scalar& relaxFactor = 0.65)
		: m_grid(grid), m_maxIteration(maxIteration), m_relaxFactor(relaxFactor) {
		m_buffer0 = std::vector<std::vector<Scalar>>(m_grid.dimension.x(), std::vector<Scalar>(m_grid.dimension.y()));
		m_buffer1 = std::vector<std::vector<Scalar>>(m_grid.dimension.x(), std::vector<Scalar>(m_grid.dimension.y()));
		m_rhs = std::vector<std::vector<Scalar>>(m_grid.dimension.x(), std::vector<Scalar>(m_grid.dimension.y()));
	}

	void BuildRHS();

	void Solve(const Scalar& liquidDensity, const Scalar& dt, bool warmStart = true);

private:
	Grid<2>& m_grid;
	std::vector<std::vector<Scalar>> m_buffer0;
	std::vector<std::vector<Scalar>> m_buffer1;
	std::vector<std::vector<Scalar>> m_rhs;
	Integer m_maxIteration;
	Scalar m_relaxFactor;

};
#pragma endregion

#pragma region Jacobi Solver 3D

#pragma endregion

#pragma endregion


#pragma region Conjugate Gradient Solver
template <int DIM>
class CGSolver;

#pragma region Conjugate Gradient Solver 2D
template <>
class CGSolver<2> {
public:
	enum class PrecondType {
		None, Jacobi, Multigrid
	};

public:
	CGSolver(Grid<2>& grid, Integer maxIteration = 500,
		Scalar tolerance = 1e-5, PrecondType type = PrecondType::Jacobi,
		Integer maxLevel = 1, Integer numPreSmooth = 5, Integer numPostSmooth = 5, Integer numFinalSmooth = 10);

	void Init();

	void BuildRHS();

	void BuildLHS(const Scalar& liquidDensity, const Scalar& dt);

	Scalar Reduce(const std::vector<std::vector<Scalar>>& vec0, const std::vector<std::vector<Scalar>>& vec1) const;

	void ApplyPreconditioner();

	Scalar GetAx(const std::vector<std::vector<Scalar>>& grid, Integer xi, Integer yi, Integer level = 0) const;

	void Solve(const Scalar& liquidDensity, const Scalar& dt);

	Scalar GetApproxLiquidPhi(Integer xi, Integer yi, Integer level) const;
	Scalar GetApproxWx(Integer xi, Integer yi, Integer level) const;
	Scalar GetApproxWy(Integer xi, Integer yi, Integer level) const;

protected:

	void Restrict(Integer toLevel);

	void Prolong(Integer toLevel);

	void Smooth(Integer smoothLevel, Integer phase);
private:
	Grid<2>& m_grid;

	std::vector<std::vector<Scalar>> m_rhs;
	std::vector<std::vector<Scalar>> m_p;
	std::vector<std::vector<Scalar>> m_x;

	std::vector<std::vector<std::vector<Scalar>>> m_r;
	std::vector<std::vector<std::vector<Scalar>>> m_z;
	std::vector<std::vector<std::vector<Scalar>>> m_Adiag;
	std::vector<std::vector<std::vector<Vector4>>> m_Acoef; // right, left, top, bottom

	Integer m_maxIteration;
	Scalar m_tolerance;
	Scalar m_oldRZ;
	Scalar m_newRZ;
	Scalar m_alpha;
	Scalar m_beta;

	PrecondType m_type;
	Integer m_maxLevel;
	Integer m_numPreSmooth;
	Integer m_numPostSmooth;
	Integer m_numFinalSmooth;

};
#pragma endregion

#pragma region Conjugate Gradient Solver 3D
template <>
class CGSolver<3> {
public:
	enum class PrecondType {
		None, Jacobi, Multigrid
	};

public:
	CGSolver(Grid<3>& grid, Integer maxIteration = 500,
		Scalar tolerance = 1e-5, PrecondType type = PrecondType::Jacobi,
		Integer maxLevel = 1, Integer numPreSmooth = 5, Integer numPostSmooth = 5, Integer numFinalSmooth = 10);

	void Init();

	void BuildRHS();

	void BuildLHS(const Scalar& liquidDensity, const Scalar& dt);

	Scalar Reduce(const std::vector<std::vector<std::vector<Scalar>>>& vec0, const std::vector<std::vector<std::vector<Scalar>>>& vec1) const;

	void ApplyPreconditioner();

	Scalar GetAx(const std::vector<std::vector<std::vector<Scalar>>>& grid, Integer xi, Integer yi, Integer zi, Integer level = 0) const;

	void Solve(const Scalar& liquidDensity, const Scalar& dt);

	Scalar GetApproxLiquidPhi(Integer xi, Integer yi, Integer zi, Integer level) const;
	Scalar GetApproxWx(Integer xi, Integer yi, Integer zi, Integer level) const;
	Scalar GetApproxWy(Integer xi, Integer yi, Integer zi, Integer level) const;
	Scalar GetApproxWz(Integer xi, Integer yi, Integer zi, Integer level) const;

protected:

	void Restrict(Integer toLevel);

	void Prolong(Integer toLevel);

	void Smooth(Integer smoothLevel, Integer phase);
private:
	Grid<3>& m_grid;

	std::vector<std::vector<std::vector<Scalar>>> m_rhs;
	std::vector<std::vector<std::vector<Scalar>>> m_p;
	std::vector<std::vector<std::vector<Scalar>>> m_x;

	std::vector<std::vector<std::vector<std::vector<Scalar>>>> m_r;
	std::vector<std::vector<std::vector<std::vector<Scalar>>>> m_z;
	std::vector<std::vector<std::vector<std::vector<Scalar>>>> m_Adiag;
	std::vector<std::vector<std::vector<std::vector<VectorX<6>>>>> m_Acoef; // right, left, top, bottom, front, back

	Integer m_maxIteration;
	Scalar m_tolerance;
	Scalar m_oldRZ;
	Scalar m_newRZ;
	Scalar m_alpha;
	Scalar m_beta;

	PrecondType m_type;
	Integer m_maxLevel;
	Integer m_numPreSmooth;
	Integer m_numPostSmooth;
	Integer m_numFinalSmooth;

};

#pragma endregion

#pragma endregion