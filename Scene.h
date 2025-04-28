#pragma once

#include "Definition.h"
#include "Boundary.h"
#include "SceneParameter.h"

template <int DIM>
class LiquidModel;

template <int DIM>
class HairModel;

template <int DIM>
class Scene {
public:
	Scene(const std::string& configFile);

	inline const SceneParameter<DIM>& GetSceneParameter() const { return m_sceneParams; }
	inline const SolverParameter<DIM>& GetSolverParameter() const { return m_solverParams; }
	inline const LiquidParameter<DIM>& GetLiquidParameter() const { return m_liquidParams; }
	inline const HairParameter<DIM>& GetHairParameter() const { return m_hairParams; }

	inline const Scalar& GetLiquidDensity() const { return m_liquidParams.density; }
	inline const Scalar& GetGravity() const { return m_sceneParams.gravity; }
	inline const VectorX<DIM>& GetMinCoord() const { return m_sceneParams.minCoord; }
	inline const VectorX<DIM>& GetMaxCoord() const { return m_sceneParams.maxCoord; }
	inline const VectorXi<DIM>& GetGridDimension() const { return m_sceneParams.gridDimension; }

	inline const Scalar& GetYoungsModulus() const { return m_hairParams.youngsModulus; }
	inline const Scalar& GetShearModulus() const { return m_hairParams.shearModulus; }

	inline LiquidModel<DIM>& GetLiquidModel() { return *m_liquidModel; }
	inline HairModel<DIM>& GetHairModel() { return *m_hairModel; }

	void PlaySimulation(const std::string& filename);

	~Scene();

private:
	SolverParameter<DIM> m_solverParams;
	SceneParameter<DIM> m_sceneParams;
	LiquidParameter<DIM> m_liquidParams;
	HairParameter<DIM> m_hairParams;

	std::unique_ptr<LiquidModel<DIM>> m_liquidModel;
	std::unique_ptr<HairModel<DIM>> m_hairModel;
	OperatorBoundary<DIM> m_boundary;

};

