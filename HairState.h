#pragma once

#include "Definition.h"
#include "Array.h"
#include "MathUtils.h"
#include "ReferenceFrame.h"
#include "Kappa.h"
#include "Twist.h"
#include "BendingProduct.h"

#include "StretchForce.h"


template <int DIM>
class HairModel;

template <int DIM>
class HairState {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	static_assert(DIM == 3, "HairState ONLY support 3D now");

	HairState(HairModel<DIM>& model, Integer hairIndex);

	void Compute(bool isInit = false);

	void ComputeBasis();

	// assume current state is the rest state, compute rest properties
	void ComputeRestProperties(const Scalar& damping = 0.0);

	void Clear();

	void AccumulateEnergy();

	void AccumulateForce();

	void AccumulateHessEnergy();

	void ApplyGravity();

	void SolveConstraints(const Scalar& dt);

	inline const VectorX<DIM>& GetPosition(Integer vertexIndex) const { return m_model.GetPosition(m_hairIndex, vertexIndex); }
	inline const VectorX<DIM>& GetVelocity(Integer vertexIndex) const { return m_model.GetVelocity(m_hairIndex, vertexIndex); }
	inline Scalar GetMass(Integer vertexIndex) const { return m_model.GetMass(m_hairIndex, vertexIndex); }
	inline Scalar GetRadius(Integer vertexIndex) const { return m_model.GetRadius(m_hairIndex, vertexIndex); }
	inline bool IsFixed(Integer vertexIndex) const { return m_model.IsFixed(m_hairIndex, vertexIndex); }

	inline void SetPosition(Integer vertexIndex, const VectorX<DIM>& position) { return m_model.SetPosition(m_hairIndex, vertexIndex, position); }
	inline void SetVelocity(Integer vertexIndex, const VectorX<DIM>& velocity) { return m_model.SetVelocity(m_hairIndex, vertexIndex, velocity); }
	inline void SetMass(Integer vertexIndex, const Scalar& mass) { return m_model.SetMass(m_hairIndex, vertexIndex, mass); }

	inline ReferenceFrame<DIM>& GetReferenceFrame() { return m_referenceFrame; }
	inline ReferenceTwist<DIM>& GetReferenceTwist() { return m_referenceTwist; }
	inline MaterialFrame<DIM>& GetMaterialFrame() { return m_materialFrame; }
	inline Kappa<DIM>& GetKappa() { return m_kappa; }
	inline GradKappa<DIM>& GetGradKappa() { return m_gradKappa; }
	inline HessKappa<DIM>& GetHessKappa() { return m_hessKappa; }
	inline Twist<DIM>& GetTwist() { return m_twist; }
	inline GradTwist<DIM>& GetGradTwist() { return m_gradTwist; }
	inline HessTwist<DIM>& GetHessTwist() { return m_hessTwist; }
	inline BendingProduct<DIM>& GetBendingProduct() { return m_bendingProduct; }
	inline Integer GetHairIndex() const { return m_hairIndex; }
	inline Integer GetNumVertices() const { return m_numVertices; }
	inline Integer GetNumEdges() const { return m_numEdges; }
	inline Integer GetNumCurvatureBinormals() const { return m_numCurvatureBinormals; }
	inline const std::vector<VectorX<DIM>>& GetKb() const { return m_kb; }
	inline const std::vector<VectorX<DIM>>& GetTangent() const { return m_tangent; }
	inline const std::vector<Scalar>& GetTheta() const { return m_theta; }
	inline const std::vector<Scalar>& GetLength() const { return m_length; }
	inline const std::vector<Matrix2>& GetBaseBendingMatrix() const { return m_baseBendingMatrix; }
	inline const HairModel<DIM>& GetModel() const { return m_model; }

	inline const Scalar& GetRestLength(Integer edgeIndex) const { return m_restLength[edgeIndex]; }
	inline const Scalar& GetVoronoiLength(Integer vertexIndex) const { return m_voronoiLength[vertexIndex]; }
	inline const Vector4& GetRestKappa(Integer edgeIndex) const { return m_restKappa[edgeIndex]; }
	inline const Scalar& GetRestTwist(Integer edgeIndex) const { return m_restTwist[edgeIndex]; }
	inline const Scalar& GetTotalRestLength() const { return m_totalRestLength; }
	inline const Scalar& GetEnergy(Integer vertexIndex) const { return m_energy[vertexIndex]; }
	inline const VectorX<DIM + 1>& GetForce(Integer vertexIndex) const { return m_force[vertexIndex]; }
	inline const VectorX<DIM>& GetKb(Integer index) const { return m_kb[index]; }
	inline const VectorX<DIM>& GetTangent(Integer index) const { return m_tangent[index]; }
	inline const Scalar& GetTheta(Integer index) const { return m_theta[index]; }
	inline const Scalar& GetLength(Integer index) const { return m_length[index]; }
	inline const Matrix2& GetBaseBendingMatrix(Integer index) const { return m_baseBendingMatrix[index]; }

private:
	HairModel<DIM>& m_model;

	Integer m_hairIndex;
	Integer m_numVertices;
	Integer m_numEdges;
	Integer m_numCurvatureBinormals;

	// basis properties
	std::vector<VectorX<DIM>> m_kb;
	std::vector<VectorX<DIM>> m_tangent;
	std::vector<Scalar> m_theta;
	std::vector<Scalar> m_length;
	std::vector<Matrix2> m_baseBendingMatrix;

	// frame, twist, and kappa properties
	ReferenceFrame<DIM> m_referenceFrame;
	ReferenceTwist<DIM> m_referenceTwist;
	MaterialFrame<DIM> m_materialFrame;
	Kappa<DIM> m_kappa;
	GradKappa<DIM> m_gradKappa;
	HessKappa<DIM> m_hessKappa;
	Twist<DIM> m_twist;
	GradTwist<DIM> m_gradTwist;
	HessTwist<DIM> m_hessTwist;
	BendingProduct<DIM> m_bendingProduct;

	// rest properties
	Scalar m_totalRestLength;
	std::vector<Scalar> m_restLength;
	std::vector<Scalar> m_voronoiLength;
	std::vector<Vector4> m_restKappa;
	std::vector<Scalar> m_restTwist;

	// force related properties
	std::vector<Scalar> m_energy;
	std::vector<VectorX<DIM + 1>> m_force; // DIM for xyz, +1 for twist theta
	std::vector<MatrixX<DIM + 1, DIM + 1>> m_hessEnergy;
	StretchForce<DIM> m_stretchForce;

};
