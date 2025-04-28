#pragma once

#include "Definition.h"
#include "MathUtils.h"


template <int DIM>
class Kappa {

public:
	Kappa(const std::vector<VectorX<DIM>>& kb, const std::vector<VectorX<DIM>>& matFrame1, const std::vector<VectorX<DIM>>& matFrame2);

	void Compute();

	void CopyFrom(const Kappa<DIM>& other);

	void Print() const;

	inline Integer Size() const { return m_size; }
	inline const Vector4& GetKappa(Integer i) const { return m_kappa[i]; }
	inline const std::vector<Vector4>& GetKappaArray() const { return m_kappa; }

private:

	const std::vector<VectorX<DIM>>& m_kb;
	const std::vector<VectorX<DIM>>& m_matFrame1;
	const std::vector<VectorX<DIM>>& m_matFrame2;
	std::vector<Vector4> m_kappa;
	Integer m_size;
};


template <int DIM>
class GradKappa {

public:
	GradKappa(const std::vector<VectorX<DIM>>& tangent, 
			  const std::vector<Scalar>& length,
			  const std::vector<VectorX<DIM>>& kb, 
		      const std::vector<VectorX<DIM>>& matFrame1, 
			  const std::vector<VectorX<DIM>>& matFrame2,
			  const std::vector<Vector4>& kappa);

	void Compute();

	void CopyFrom(const GradKappa<DIM>& other);

	void Print() const;

	inline Integer Size() const { return m_size; }
	inline const MatrixX<11, 4>& GetGradKappa(Integer i) const { return m_gradKappa[i]; }
	inline const std::vector<MatrixX<11, 4>>& GetGradKappaArray() const { return m_gradKappa; }

private:

	const std::vector<VectorX<DIM>>& m_tangent;
	const std::vector<Scalar>& m_length;
	const std::vector<VectorX<DIM>>& m_kb;
	const std::vector<VectorX<DIM>>& m_matFrame1;
	const std::vector<VectorX<DIM>>& m_matFrame2;
	const std::vector<Vector4>& m_kappa;
	std::vector<MatrixX<11, 4>> m_gradKappa;
	Integer m_size;
};

template <int DIM>
class HessKappa {
public:
	HessKappa(const std::vector<VectorX<DIM>>& tangent,
		const std::vector<Scalar>& length,
		const std::vector<VectorX<DIM>>& kb,
		const std::vector<VectorX<DIM>>& matFrame1,
		const std::vector<VectorX<DIM>>& matFrame2,
		const std::vector<Vector4>& kappa);

	void Compute();

	void CopyFrom(const HessKappa<DIM>& other);

	void Print() const;

	inline Integer Size() const { return m_size; }
	inline const MatrixX<11, 11>& GetHessKappa(Integer i) const { return m_hessKappa[i]; }
	inline const std::vector<MatrixX<11, 11>>& GetHessKappaArray() const { return m_hessKappa; }

private:
	const std::vector<VectorX<DIM>>& m_tangent;
	const std::vector<Scalar>& m_length;
	const std::vector<VectorX<DIM>>& m_kb;
	const std::vector<VectorX<DIM>>& m_matFrame1;
	const std::vector<VectorX<DIM>>& m_matFrame2;
	const std::vector<Vector4>& m_kappa;
	std::vector<MatrixX<11, 11>> m_hessKappa;
	Integer m_size;

};