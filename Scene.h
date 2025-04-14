#pragma once

#include "Definition.h"
#include "Boundary.h"
#include "SceneParameter.h"

template <int DIM>
class LiquidModel;

template <int DIM>
class Scene {
public:
	Scene(SceneParameter<DIM> params);

	inline const SceneParameter<DIM>& GetParameter() const { return m_params; }
	inline const Scalar& GetLiquidDensity() const { return m_params.liquidDensity; }
	inline const Scalar& GetGravity() const { return m_params.gravity; }
	inline const Integer& GetNumParticle() const { return m_params.numParticles; }
	inline const VectorX<DIM>& GetMinCoord() const { return m_params.minCoord; }
	inline const VectorX<DIM>& GetMaxCoord() const { return m_params.maxCoord; }
	inline const VectorXi<DIM>& GetGridDimension() const { return m_params.gridDimension; }

	inline LiquidModel<DIM>& GetLiquidModel() { return m_model; }


private:
	SceneParameter<DIM> m_params;
	LiquidModel<DIM> m_model;
	OperatorBoundary<DIM> m_boundary;

};

