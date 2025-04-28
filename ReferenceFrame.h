#pragma once

#include "Definition.h"
#include "MathUtils.h"

template <int DIM>
class ReferenceFrame {
public:
    ReferenceFrame(const std::vector<VectorX<DIM>>& tangent)
        : m_tangent(tangent), m_prevTangent(tangent), m_size(static_cast<Integer>(tangent.size())) {
        m_refFrame1.resize(m_size);
        m_refFrame2.resize(m_size);
    }

    void Init() {
        m_refFrame1[0] = MathUtils::ComputeNormal<DIM>(m_tangent[0]);
        m_refFrame2[0] = (m_tangent[0].cross(m_refFrame1[0])).normalized();
        for (Integer i = 1; i < m_size; ++i) {
            m_refFrame1[i] = MathUtils::ParallelTransport(m_refFrame1[i - 1], m_tangent[i - 1], m_tangent[i]);
            m_refFrame2[i] = MathUtils::ParallelTransport(m_refFrame2[i - 1], m_tangent[i - 1], m_tangent[i]);
        }
    }

    void Compute() {
#pragma omp parallel for
        for (Integer i = 0; i < m_size; ++i) {
            VectorX<DIM> ref1 = MathUtils::ParallelTransport(
                m_refFrame1[i], m_prevTangent[i], m_tangent[i]);
            m_refFrame1[i] = MathUtils::OrthoNormalized(ref1, m_tangent[i]);
            m_refFrame2[i] = (m_tangent[i].cross(m_refFrame1[i])).normalized();
            m_prevTangent[i] = m_tangent[i];
        }
    }

    void CopyFrom(const ReferenceFrame<DIM>& other) {
#ifdef _DEBUG
        ASSERT_MSG(m_size == other.m_size, "Size mismatch!");
#endif
        for (Integer i = 0; i < m_size; ++i) {
            m_refFrame1[i] = other.m_refFrame1[i];
            m_refFrame2[i] = other.m_refFrame2[i];
            m_prevTangent[i] = other.m_prevTangent[i];
        }
    }

    void Print() const {
        for (Integer i = 0; i < m_size; ++i) {
            std::cout << "RefFrame1[" << i << "]: "
                << m_refFrame1[i].transpose() << "\n";
            std::cout << "RefFrame2[" << i << "]: "
                << m_refFrame2[i].transpose() << "\n";
        }
    }

    inline Integer Size() const { return m_size; }
    inline const VectorX<DIM>& GetRefFrame1(Integer i) const { return m_refFrame1[i]; }
    inline const VectorX<DIM>& GetRefFrame2(Integer i) const { return m_refFrame2[i]; }
    inline const VectorX<DIM>& GetPrevTangent(Integer i) const { return m_prevTangent[i]; }

    inline const std::vector<VectorX<DIM>>& GetRefFrame1Array() const { return m_refFrame1; }
    inline const std::vector<VectorX<DIM>>& GetRefFrame2Array() const { return m_refFrame2; }
    inline const std::vector<VectorX<DIM>>& GetPrevTangentArray() const { return m_prevTangent; }

private:
    Integer m_size;
    const std::vector<VectorX<DIM>>& m_tangent;
    std::vector<VectorX<DIM>> m_refFrame1;
    std::vector<VectorX<DIM>> m_refFrame2;
    std::vector<VectorX<DIM>> m_prevTangent;
};


template <int DIM>
class MaterialFrame {
public:
	MaterialFrame(const std::vector<Scalar>& theta, const std::vector<VectorX<DIM>>& refFrame1, const std::vector<VectorX<DIM>>& refFrame2)
		: m_refFrame1(refFrame1), m_refFrame2(refFrame2), m_theta(theta), m_size(static_cast<Integer>(refFrame1.size())) {
        m_materialFrame1.resize(m_size);
		m_materialFrame2.resize(m_size);
	}

    void Compute() {
#pragma omp parallel for
        for (Integer index = 0; index < m_size; ++index) {
			Scalar ct = std::cos(m_theta[index]);
            Scalar st = std::sin(m_theta[index]);
			m_materialFrame1[index] = ct * m_refFrame1[index] + st * m_refFrame2[index];
			m_materialFrame2[index] = -st * m_refFrame1[index] + ct * m_refFrame2[index];
        }
    }

    void CopyFrom(const MaterialFrame<DIM>& other) {
#ifdef _DEBUG
        ASSERT_MSG(m_size == other.m_size, "Size mismatch!");
#endif

		for (Integer i = 0; i < m_size; ++i) {
			m_materialFrame1[i] = other.m_materialFrame1[i];
			m_materialFrame2[i] = other.m_materialFrame2[i];
		}
    }

	void Print() const {
		for (Integer i = 0; i < m_size; ++i) {
			std::cout << "MaterialFrame1[" << i << "]: " << m_materialFrame1[i].transpose() << "\n";
			std::cout << "MaterialFrame2[" << i << "]: " << m_materialFrame2[i].transpose() << "\n";
		}
	}

	inline Integer Size() const { return m_size; }
	inline const VectorX<DIM>& GetMaterialFrame1(Integer i) const { return m_materialFrame1[i]; }
	inline const VectorX<DIM>& GetMaterialFrame2(Integer i) const { return m_materialFrame2[i]; }
	inline const std::vector<VectorX<DIM>>& GetMaterialFrame1Array() const { return m_materialFrame1; }
	inline const std::vector<VectorX<DIM>>& GetMaterialFrame2Array() const { return m_materialFrame2; }

private:
	const std::vector<VectorX<DIM>>& m_refFrame1;
	const std::vector<VectorX<DIM>>& m_refFrame2;
	const std::vector<Scalar>& m_theta;
	std::vector<VectorX<DIM>> m_materialFrame1;
    std::vector<VectorX<DIM>> m_materialFrame2;
	Integer m_size;
};


template <int DIM>
class ReferenceTwist {
public:
    ReferenceTwist(const std::vector<VectorX<DIM>>& tangent, const std::vector<VectorX<DIM>>& refFrame1, const std::vector<VectorX<DIM>>& refFrame2)
        : m_tangent(tangent), m_refFrame1(refFrame1), m_refFrame2(refFrame2), m_size(static_cast<Integer>(refFrame1.size()) - 1) {
        m_refTwist.resize(m_size);
    }

    void Compute() {
#pragma omp parallel for
        for (Integer i = 0; i < m_size; ++i) {
            const VectorX<DIM>& prevF = m_refFrame1[i];
            const VectorX<DIM>& nextF = m_refFrame1[i + 1];
            VectorX<DIM> pred = MathUtils::ParallelTransport(
                prevF, m_tangent[i], m_tangent[i + 1]);
            pred = MathUtils::Rotate(pred, m_tangent[i + 1], m_refTwist[i]);
            m_refTwist[i] = MathUtils::ComputeAngle(pred, nextF, m_tangent[i + 1]);
        }
    }

    void CopyFrom(const ReferenceTwist<DIM>& other) {
#ifdef _DEBUG
		ASSERT_MSG(m_size == other.m_size, "Size mismatch!");
#endif
        for (Integer i = 0; i < m_size; ++i) {
            m_refTwist[i] = other.m_refTwist[i];
        }
    }

    void Print() const {
        for (Integer i = 0; i < m_size; ++i) {
            std::cout << "RefTwist[" << i << "]: " << m_refTwist[i] << "\n";
        }
    }

    inline Integer Size() const { return m_size; }
	inline const Scalar& GetRefTwist(Integer i) const { return m_refTwist[i]; }
    inline const std::vector<Scalar>& GetRefTwistArray() const { return m_refTwist; }

private:
    const std::vector<VectorX<DIM>>& m_tangent;
    const std::vector<VectorX<DIM>>& m_refFrame1;
    const std::vector<VectorX<DIM>>& m_refFrame2;
    std::vector<Scalar> m_refTwist;
    Integer m_size;
};