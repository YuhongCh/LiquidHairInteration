#pragma once

#include "Definition.h"
#include "MathUtils.h"

template <int DIM>
class Twist {
public:
    Twist(const std::vector<Scalar>& theta, const std::vector<Scalar>& refTwist)
        : m_theta(theta), m_refTwist(refTwist), m_size(static_cast<Integer>(refTwist.size())) {
        m_twist.resize(m_size);
    }

    void Compute() {
#pragma omp parallel for
		for (Integer i = 0; i < m_size; ++i) {
			const Scalar& prevTheta = m_theta[i];
			const Scalar& nextTheta = m_theta[i + 1];
			m_twist[i] = m_refTwist[i] + nextTheta - prevTheta;
		}
    }

	void CopyFrom(const Twist& other) {
#ifdef _DEBUG
		ASSERT_MSG(m_size == other.m_size, "Size mismatch!");
#endif
		for (Integer i = 0; i < m_size; ++i) {
			m_twist[i] = other.m_twist[i];
		}
	}
	void Print() const {
		for (Integer i = 0; i < m_size; ++i) {
			std::cout << "Twist[" << i << "]: " << m_twist[i] << "\n";
		}
	}
	inline Integer Size() const { return m_size; }
	inline const Scalar& GetTwist(Integer i) const { return m_twist[i]; }
	inline const std::vector<Scalar>& GetTwistArray() const { return m_twist; }


private:
    const std::vector<Scalar>& m_theta;
    const std::vector<Scalar>& m_refTwist;
    std::vector<Scalar> m_twist;
    Integer m_size;
};

template <int DIM>
class GradTwist {
public:
	GradTwist(const std::vector<Scalar>& length, const std::vector<VectorX<DIM>>& kb)
		: m_length(length), m_kb(kb), m_size(static_cast<Integer>(kb.size())) {
		m_gradTwist.resize(m_size);
	}

	void Compute() {
#pragma omp parallel for
		for (Integer index = 0; index < m_size; ++index) {
			const VectorX<DIM>& kb = m_kb[index];
			const Scalar& prevLength = m_length[index];
			const Scalar& nextLength = m_length[index + 1];
			VectorX<11>& gradTwist = m_gradTwist[index];

			gradTwist.segment<3>(0) = -0.5 / prevLength * kb;
			gradTwist.segment<3>(8) = 0.5 / nextLength * kb;
			gradTwist.segment<3>(4) = -(gradTwist.segment<3>(0) + gradTwist.segment<3>(8));
			gradTwist(3) = -1;
			gradTwist(7) = 1;
		}
	}

	void CopyFrom(const GradTwist& other) {
#ifdef _DEBUG
		ASSERT_MSG(m_size == other.m_size, "Size mismatch!");
#endif
		for (Integer i = 0; i < m_size; ++i) {
			m_gradTwist[i] = other.m_gradTwist[i];
		}
	}
	void Print() const {
		for (Integer i = 0; i < m_size; ++i) {
			std::cout << "GradTwist[" << i << "]: " << m_gradTwist[i].transpose() << "\n";
		}
	}

	inline Integer Size() const { return m_size; }
	inline const VectorX<11>& GetGradTwist(Integer i) const { return m_gradTwist[i]; }
	inline const std::vector<VectorX<11>>& GetGradTwistArray() const { return m_gradTwist; }

private:
	const std::vector<Scalar>& m_length;
	const std::vector<VectorX<DIM>>& m_kb;
	std::vector<VectorX<11>> m_gradTwist;
	Integer m_size;
};


template <int DIM>
class HessTwist {
public:
	HessTwist(const std::vector<VectorX<DIM>>& tangent, const std::vector<Scalar>& length, const std::vector<VectorX<DIM>>& kb)
		: m_tangent(tangent), m_length(length), m_kb(kb), m_size(static_cast<Integer>(kb.size())) {
		m_hessTwist.resize(m_size);
	}

	void Compute() {
#pragma omp parallel for
		for (Integer index = 0; index < m_size; ++index) {
			const Scalar& ls = m_length[index];
			const Scalar& le = m_length[index + 1];
			const VectorX<DIM>& ts = m_tangent[index];
			const VectorX<DIM>& te = m_tangent[index + 1];
			const VectorX<DIM>& kb = m_kb[index];
			MatrixX<11, 11>& hessTwist = m_hessTwist[index];

			Scalar chi = MathUtils::Clamp(1.0 + ts.dot(te), EPSILON, INF);
			Scalar inv_chi = 1.0 / chi;
			VectorX<DIM> tilde_t = inv_chi * (ts + te);
			MatrixX<DIM, DIM> D2mDs2 = -0.25 / (ls * ls) * (kb * (ts + tilde_t).transpose() + (ts + tilde_t) * kb.transpose());
			MatrixX<DIM, DIM> D2mDe2 = -0.25 / (le * le) * (kb * (te + tilde_t).transpose() + (te + tilde_t) * kb.transpose());
			MatrixX<DIM, DIM> D2mDsDe = 0.5 / (ls * le) * (2.0 / chi * MathUtils::CrossMatrix(ts) - kb * tilde_t.transpose());
			MatrixX<DIM, DIM> D2mDeDs = D2mDsDe.transpose();

			hessTwist.block<3, 3>(0, 0) = D2mDs2;
			hessTwist.block<3, 3>(0, 4) = -D2mDs2 + D2mDsDe;
			hessTwist.block<3, 3>(4, 0) = -D2mDs2 + D2mDeDs;
			hessTwist.block<3, 3>(4, 4) = D2mDs2 - (D2mDsDe + D2mDeDs) + D2mDe2;
			hessTwist.block<3, 3>(0, 8) = -D2mDsDe;
			hessTwist.block<3, 3>(8, 0) = -D2mDeDs;
			hessTwist.block<3, 3>(4, 8) = D2mDsDe - D2mDe2;
			hessTwist.block<3, 3>(8, 4) = D2mDeDs - D2mDe2;
			hessTwist.block<3, 3>(8, 8) = D2mDe2;
		}
	}

	void CopyFrom(const HessTwist& other) {
#ifdef _DEBUG
		ASSERT_MSG(m_size == other.m_size, "Size mismatch!");
#endif
#pragma omp simd
		for (Integer i = 0; i < m_size; ++i) {
			m_hessTwist[i] = other.m_hessTwist[i];
		}
	}
	void Print() const {
		for (Integer i = 0; i < m_size; ++i) {
			std::cout << "HessTwist[" << i << "]: " << m_hessTwist[i].transpose() << "\n";
		}
	}
	inline Integer Size() const { return m_size; }
	inline const MatrixX<11, 11>& GetHessTwist(Integer i) const { return m_hessTwist[i]; }
	inline const std::vector<MatrixX<11, 11>>& GetHessTwistArray() const { return m_hessTwist; }

private:
	const std::vector<VectorX<DIM>>& m_tangent;
	const std::vector<Scalar>& m_length;
	const std::vector<VectorX<DIM>>& m_kb;
	std::vector<MatrixX<11, 11>> m_hessTwist;
	Integer m_size;

};
