#pragma once

#include "Grid.h"
#include "Array.h"

template <int DIM>
class SolverParameter;

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
		Integer maxLevel = 1, Integer numPreSmooth = 5, Integer numPostSmooth = 5, Integer numFinalSmooth = 10, Scalar smoothFactor = 0.5);

	CGSolver(Grid<2>& grid, const SolverParameter<2>& params);

	void Init();

	void BuildRHS();

	void BuildLHS(const Scalar& liquidDensity, const Scalar& dt);

	void ApplyPreconditioner();

	Scalar GetAx(const Array2<Scalar>& grid, Integer xi, Integer yi, Integer level = 0) const;

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

	Array2<Scalar> m_rhs;
	Array2<Scalar> m_p;
	Array2<Scalar> m_x;

	std::vector<Array2<Scalar>> m_r;
	std::vector<Array2<Scalar>> m_z;
	std::vector<Array2<Scalar>> m_Adiag;
	std::vector<Array2<Vector4>> m_Acoef; // right, left, top, bottom

	Integer m_maxIteration;
	Scalar m_smoothFactor;
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
		Integer maxLevel = 1, Integer numPreSmooth = 5, Integer numPostSmooth = 5, Integer numFinalSmooth = 10, Scalar smoothFactor = 0.5);

	CGSolver(Grid<3>& grid, const SolverParameter<3>& params);

	void Init();

	void BuildRHS();

	void BuildLHS(const Scalar& liquidDensity, const Scalar& dt);

	void ApplyPreconditioner();

	Scalar GetAx(const Array3<Scalar>& grid, Integer xi, Integer yi, Integer zi, Integer level = 0) const;

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

	 Array3<Scalar> m_rhs;
	 Array3<Scalar> m_p;
	 Array3<Scalar> m_x;

	std::vector<Array3<Scalar>> m_r;
	std::vector<Array3<Scalar>> m_z;
	std::vector<Array3<Scalar>> m_Adiag;
	std::vector<Array3<VectorX<6>>> m_Acoef; // right, left, top, bottom, front, back

	Integer m_maxIteration;
	Scalar m_smoothFactor;
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