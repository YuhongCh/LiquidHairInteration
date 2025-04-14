#pragma once

#include "Definition.h"
#include "Grid.h"
#include "ParticleSystem.h"
#include "Boundary.h"
#include "Scene.h"
#include "NeighborSearch.h"
#include "PressureSolver.h"

template <int DIM>
class LiquidModel {
public:
	LiquidModel(const Scene<DIM>& scene);

	void Init();

	void PreStep(const Scalar& dt);

	void Step(const Scalar& dt);

	void PostStep(const Scalar& dt);

public:
	void ComputeLiquidPhi();

	void ComputeSolidPhi();

	void Particle2Grid();

	void Grid2Particle();

	void Projection(const Scalar& dt);

	void SolvePressure(const Scalar& dt);

	void ApplyGravity(const Scalar& dt);

	void AdvectParticle(const Scalar& dt);

	void Extrapolate();

	void CorrectVolume(Scalar dt);

	Scalar ComputeDivergence() const;

	inline Scalar GetCFL() const {
		Scalar maxVel = 0.0;
		for (Integer pi = 0; pi < m_ps.NumParticles(); ++pi) {
			maxVel = max(maxVel, m_ps.GetVelocity(pi).cwiseAbs().maxCoeff());
		}
		return MathUtils::Clamp(m_grid.dx.x() / maxVel, 0.001, 0.1);
	}

	inline void SetBoundary(BaseBoundary<DIM>* boundary) { m_boundary = boundary; }

	inline const ParticleSystem<DIM>& GetParticleSystem() const { return m_ps; }

private:

	Grid<DIM> m_grid;
	ParticleSystem<DIM> m_ps;
	HashNeighborSearch<DIM> m_ns;
	BaseBoundary<DIM>* m_boundary;
	const Scene<DIM>& m_scene;
	CGSolver<DIM> m_solver;

	Integer m_correctStep;
	Integer m_correctCycle;
};

