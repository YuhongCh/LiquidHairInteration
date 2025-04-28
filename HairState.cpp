#include "HairState.h"
#include "HairModel.h"

//template class HairState<3>;

template <>
HairState<3>::HairState(HairModel<3>& model, Integer hairIndex)
	: m_model(model), m_hairIndex(hairIndex), m_numVertices(model.NumVertexOnHair(hairIndex)),
	m_numEdges(m_numVertices - 1), m_numCurvatureBinormals(m_numEdges - 1),
	m_kb(m_numCurvatureBinormals), m_tangent(m_numEdges), m_theta(m_numVertices),
	m_length(m_numEdges), m_baseBendingMatrix(m_numVertices), m_referenceFrame(m_tangent),
	m_referenceTwist(m_tangent, m_referenceFrame.GetRefFrame1Array(), m_referenceFrame.GetRefFrame2Array()),
	m_materialFrame(m_theta, m_referenceFrame.GetRefFrame1Array(), m_referenceFrame.GetRefFrame2Array()),
	m_kappa(m_kb, m_materialFrame.GetMaterialFrame1Array(), m_materialFrame.GetMaterialFrame2Array()),
	m_gradKappa(m_tangent, m_length, m_kb, m_materialFrame.GetMaterialFrame1Array(), m_materialFrame.GetMaterialFrame2Array(), m_kappa.GetKappaArray()),
	m_hessKappa(m_tangent, m_length, m_kb, m_materialFrame.GetMaterialFrame1Array(), m_materialFrame.GetMaterialFrame2Array(), m_kappa.GetKappaArray()),
	m_twist(m_theta, m_referenceTwist.GetRefTwistArray()), 
	m_gradTwist(m_length, m_kb),
	m_hessTwist(m_tangent, m_length, m_kb),
	m_bendingProduct(m_baseBendingMatrix, m_gradKappa.GetGradKappaArray()),
	m_restLength(m_numEdges), m_voronoiLength(m_numVertices), m_restKappa(m_numEdges), m_restTwist(m_numEdges),
	m_energy(m_numVertices), m_force(m_numVertices),
	m_stretchForce(*this) {

	// Compute baseBendingMatrix
	for (Integer vertexIndex = 0; vertexIndex < m_numVertices; ++vertexIndex) {
		const Scalar& radius = m_model.GetRadius(m_hairIndex, vertexIndex);
		m_baseBendingMatrix[vertexIndex] = std::pow(radius, 4) * PI_4 * Matrix2::Identity();
	}
	ComputeRestProperties(0.0);
	Compute();
}

template <>
void HairState<3>::ComputeRestProperties(const Scalar& damping) {
	m_totalRestLength = 0.0;
	for (Integer edgeIndex = 0; edgeIndex < m_numEdges; ++edgeIndex) {
		m_restLength[edgeIndex] = MathUtils::Lerp(m_length[edgeIndex], m_restLength[edgeIndex], damping);
		m_totalRestLength += m_restLength[edgeIndex];
	}

	m_voronoiLength[0] = 0.5 * m_restLength[0];
	m_voronoiLength[m_numEdges] = 0.5 * m_restLength[m_numEdges - 1];
	for (Integer vertexIndex = 1; vertexIndex < m_numEdges; ++vertexIndex) {
		m_voronoiLength[vertexIndex] = 0.5 * (m_restLength[vertexIndex - 1] + m_restLength[vertexIndex]);
	}

	for (Integer vertexIndex = 0; vertexIndex < m_numVertices; ++vertexIndex) {
		const Scalar& radius = GetRadius(vertexIndex);
		SetMass(vertexIndex, PI * radius * radius * m_voronoiLength[vertexIndex]);
	}

	for (Integer kbIndex = 0; kbIndex < m_numCurvatureBinormals; ++kbIndex) {
		m_restTwist[kbIndex] = MathUtils::Lerp(m_twist.GetTwist(kbIndex), m_restTwist[kbIndex], damping);
		m_restKappa[kbIndex] = MathUtils::Lerp(m_kappa.GetKappa(kbIndex), m_restKappa[kbIndex], damping);
	}
}

template <>
void HairState<3>::ComputeBasis() {
#pragma omp parallel for
	for (Integer edgeIndex = 0; edgeIndex < m_numEdges; ++edgeIndex) {
		const Vector3& p0 = m_model.GetPosition(m_hairIndex, edgeIndex);
		const Vector3& p1 = m_model.GetPosition(m_hairIndex, edgeIndex + 1);
		m_tangent[edgeIndex] = p1 - p0;
		m_length[edgeIndex] = m_tangent[edgeIndex].norm();
		m_tangent[edgeIndex] /= m_length[edgeIndex];
	}

#pragma omp parallel for
	for (Integer kbIndex = 0; kbIndex < m_numCurvatureBinormals; ++kbIndex) {
		const Vector3& ts = m_tangent[kbIndex];
		const Vector3& te = m_tangent[kbIndex + 1];
		Scalar denominator = 1.0 + te.dot(ts);
		if (MathUtils::IsSmall(denominator)) {
			m_kb[kbIndex] = 4.0 * std::tan(0.5 * std::acos(denominator - 1.0)) * MathUtils::ComputeNormal(te);
		}
		else {
			m_kb[kbIndex] = 2 * ts.cross(te) / denominator;
		}
	}
}

template <>
void HairState<3>::Compute(bool isInit) {
	ComputeBasis();
	if (isInit) m_referenceFrame.Init();
	else m_referenceFrame.Compute();

	m_referenceTwist.Compute();
	m_materialFrame.Compute();
	m_kappa.Compute();
	m_gradKappa.Compute();
	//m_hessKappa.Compute();
	m_twist.Compute();
	m_gradTwist.Compute();
	m_hessTwist.Compute();
	m_bendingProduct.Compute();
}

template <>
void HairState<3>::Clear() {
#pragma omp simd
	for (Integer vertexIndex = 0; vertexIndex < m_numVertices; ++vertexIndex) {
		m_energy[vertexIndex] = 0.0;
		m_force[vertexIndex].setZero();
		m_hessEnergy[vertexIndex].setZero();
	}
}

template <>
void HairState<3>::ApplyGravity() {
	Vector3 gravity(0.0, m_model.GetScene().GetGravity(), 0.0);

#pragma omp parallel for
	for (Integer vertexIndex = 0; vertexIndex < m_numVertices; ++vertexIndex) {
		// accumulate energy
		m_energy[vertexIndex] += gravity.dot(GetPosition(vertexIndex));

		// accumulate force
		m_force[vertexIndex].segment<3>(0) += gravity * GetMass(vertexIndex);

		// dont need accumulate hessian for gravity
	}

}

template <>
void HairState<3>::AccumulateEnergy() {
	m_stretchForce.AccumulateEnergy(m_energy);
}

template <>
void HairState<3>::AccumulateForce() {
	m_stretchForce.AccumulateForce(m_force);
}

template <>
void HairState<3>::AccumulateHessEnergy() {
	m_stretchForce.AccumulateHessEnergy(m_hessEnergy);
}