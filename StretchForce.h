#pragma once

#include "Definition.h"

template <int DIM>
class HairState;


template <int DIM>
class StretchForce {
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
		static_assert(DIM == 3, "StretchForce ONLY support 3D now");
public:
	StretchForce(const HairState<DIM>& hairState);

	void ComputeEnergyLocal(Integer index, Scalar&) const;

	void ComputeForceLocal(Integer index, VectorX<DIM>&) const;

	void ComputeForceJacobiLocal(Integer index, MatrixX<DIM, DIM>&) const;

	void AccumulateEnergy(std::vector<Scalar>&) const;

	void AccumulateForce(std::vector<VectorX<DIM + 1>>& force) const;

	void AccumulateHessEnergy(std::vector<MatrixX<DIM + 1, DIM + 1>>& jacobi) const;

private:
	const HairState<DIM>& m_hairState;

};

