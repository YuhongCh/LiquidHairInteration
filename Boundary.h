#pragma once

#include "Definition.h"

template <int DIM>
class BaseBoundary {
public:
	BaseBoundary()
		: m_leftChild(nullptr), m_rightChild(nullptr){}

	virtual Scalar Compute(const VectorX<DIM>& pos) const = 0;

	inline void SetChild(BaseBoundary* boundary) { m_leftChild = boundary; }
	inline void SetLeftChild(BaseBoundary* boundary) { m_leftChild = boundary; }
	inline void SetRightChild(BaseBoundary* boundary) { m_rightChild = boundary; }

	~BaseBoundary() {
		if (this->m_leftChild != nullptr) delete m_leftChild;
		if (this->m_rightChild != nullptr) delete m_rightChild;
	}

protected:
	BaseBoundary* m_leftChild;
	BaseBoundary* m_rightChild;
};

template <int DIM>
class OperatorBoundary : public BaseBoundary<DIM> {
public:
	enum class OperatorType {
		Union, Intersection, Substract, None
	};

public:
	OperatorBoundary(OperatorType type = OperatorType::None)
		: BaseBoundary<DIM>(), m_type(type) {}

	Scalar Compute(const VectorX<DIM>& pos) const override {
		switch (m_type) {
			case OperatorType::Union: {
				return min(this->m_leftChild->Compute(pos), this->m_rightChild->Compute(pos));
			}
			case OperatorType::Intersection: {
				return max(this->m_leftChild->Compute(pos), this->m_rightChild->Compute(pos));
			}
			case OperatorType::Substract: {
				return max(this->m_leftChild->Compute(pos), -this->m_rightChild->Compute(pos));
			}
			case OperatorType::None: {
				return this->m_leftChild->Compute(pos);
			}
			default: {
				throw std::runtime_error("Failed to recognize OperatorType");
			}
		}
		return 1.0;
	}

protected:
	OperatorType m_type;

};

template <int DIM>
class SolidBoundary : public BaseBoundary<DIM> {
public:
	SolidBoundary(const VectorX<DIM>& center, bool useInside = false)
		: BaseBoundary<DIM>(), m_center(center), m_useInside(useInside) {}

	Scalar Compute(const VectorX<DIM>& pos) const override { return 1.0; }

protected:
	VectorX<DIM> m_center;
	bool m_useInside;
};

template <int DIM>
class RectBoundary : public SolidBoundary<DIM> {
public:
	RectBoundary(const VectorX<DIM>& center, const VectorX<DIM>& halfExtent, bool useInside = false)
		: SolidBoundary<DIM>(center, useInside), m_halfExtent(halfExtent) {}

	Scalar Compute(const VectorX<DIM>& pos) const override {
		VectorX<DIM> p = (pos - this->m_center).cwiseAbs() - m_halfExtent;
		VectorX<DIM> pmax = p.cwiseMax(Scalar(0));
		Scalar outsideDist = pmax.norm();
		Scalar insideDist = min(p.maxCoeff(), Scalar(0));

		return (this->m_useInside) ? -outsideDist - insideDist : outsideDist + insideDist;
	}

protected:
	VectorX<DIM> m_halfExtent;
};


