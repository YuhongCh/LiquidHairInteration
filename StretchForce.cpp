#include "StretchForce.h"
#include "HairState.h"
#include "HairModel.h"

template class StretchForce<3>;

template <>
StretchForce<3>::StretchForce(const HairState<3>& hairState)
	: m_hairState(hairState) {
}

template <>
void StretchForce<3>::ComputeEnergyLocal(Integer index, Scalar& energy) const {
	Scalar length = m_hairState.GetLength(index);
	Scalar refLength = m_hairState.GetRestLength(index);
	Scalar youngsModulus = m_hairState.GetModel().GetScene().GetYoungsModulus();
	energy = length / refLength - 1.0;
	energy = 0.5 * youngsModulus * energy * energy * refLength;
}

template <>
void StretchForce<3>::ComputeForceLocal(Integer index, VectorX<3>& force) const {
	Scalar youngsModulus = m_hairState.GetModel().GetScene().GetYoungsModulus();
	Scalar length = m_hairState.GetLength(index);
	Scalar refLength = m_hairState.GetRestLength(index);
	Vector3 tangent = m_hairState.GetTangent(index);
	force = youngsModulus * (length / refLength - 1.0) * tangent;
}

template <>
void StretchForce<3>::ComputeForceJacobiLocal(Integer index, MatrixX<3, 3>& jacobi) const {
	Scalar youngsModulus = m_hairState.GetModel().GetScene().GetYoungsModulus();
	Scalar length = m_hairState.GetLength(index);
	Scalar refLength = m_hairState.GetRestLength(index);
	Vector3 tangent = m_hairState.GetTangent(index);
	jacobi = youngsModulus * (1.0 / refLength) * tangent * tangent.transpose();
}

template <>
void StretchForce<3>::AccumulateEnergy(std::vector<Scalar>& energy) const {
#pragma omp parallel for
	for (Integer index = 0; index < m_hairState.GetNumEdges(); ++index) {
		Scalar localEnergy;
		ComputeEnergyLocal(index, localEnergy);
		energy[index] += localEnergy;
	}
}

template <>
void StretchForce<3>::AccumulateForce(std::vector<VectorX<4>>& force) const {
#pragma omp parallel for
	for (Integer index = 0; index < m_hairState.GetNumEdges(); ++index) {
		Vector3 localForce;
		ComputeForceLocal(index, localForce);
		force[index].segment<3>(0) += localForce;
		force[index + 1].segment<3>(0) -= localForce;
	}
}

template <>
void StretchForce<3>::AccumulateHessEnergy(std::vector<MatrixX<4, 4>>& jacobi) const {
#pragma omp parallel for
	for (Integer index = 0; index < m_hairState.GetNumEdges(); ++index) {
		Matrix3 localJacobi;
		ComputeForceJacobiLocal(index, localJacobi);
		jacobi[index].block<3, 3>(0, 0) += localJacobi;
	}
}

