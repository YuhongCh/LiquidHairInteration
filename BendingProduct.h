#pragma once

#include "Definition.h"
#include "MathUtils.h"

template <int DIM>
class BendingProduct {
public:
    BendingProduct(const std::vector<Matrix2>& baseBendingMatrix, const std::vector<MatrixX<11, 4>>& gradKappa)
        : m_gradKappa(gradKappa), m_baseBendingMatrix(baseBendingMatrix), m_size(static_cast<Integer>(gradKappa.size())) {
		m_bendingProduct.resize(m_size);
    }

    void Compute() {
//#pragma omp parallel for
//        for (Integer index = 0; index < m_size; ++index) {
//            const Matrix2& bendingMat = m_baseBendingMatrix[index];
//			const MatrixX<11, 4>& gradKappa = m_gradKappa[index];
//			MatrixX<11, 11>& bendingProduct = m_bendingProduct[index];
//
//			bendingProduct = gradKappa.block<11, 2>(0, 0) * bendingMat * gradKappa.block<11, 2>(0, 0).transpose() +
//							 gradKappa.block<11, 2>(0, 2) * bendingMat * gradKappa.block<11, 2>(0, 2).transpose();
//        }
    }

	void CopyFrom(const BendingProduct& other) {
#ifdef _DEBUG
		ASSERT_MSG(m_size == other.m_size, "Size mismatch!");
#endif
		for (Integer i = 0; i < m_size; ++i) {
			m_bendingProduct[i] = other.m_bendingProduct[i];
		}
	}
	void Print() const {
		for (Integer i = 0; i < m_size; ++i) {
			std::cout << "BendingProduct[" << i << "]: " << m_bendingProduct[i].transpose() << "\n";
		}
	}
	inline Integer Size() const { return m_size; }
	inline const MatrixX<11, 11>& GetBendingProduct(Integer i) const { return m_bendingProduct[i]; }
	inline const std::vector<MatrixX<11, 11>>& GetBendingProductArray() const { return m_bendingProduct; }


private:
	const std::vector<MatrixX<11, 4>>& m_gradKappa;
	const std::vector<Matrix2>& m_baseBendingMatrix;
	std::vector<MatrixX<11, 11>> m_bendingProduct;
	Integer m_size;
};

