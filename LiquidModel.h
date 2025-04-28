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

	~LiquidModel();

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

	inline ParticleSystem<DIM>& GetParticleSystem() { return m_ps; }

	void Render() const;

private:

	Grid<DIM> m_grid;
	ParticleSystem<DIM> m_ps;
	HashNeighborSearch<DIM> m_ns;
	BaseBoundary<DIM>* m_boundary;
	const Scene<DIM>& m_scene;
	CGSolver<DIM> m_solver;

	Integer m_correctCycle;
	Integer m_correctStep;

	GLuint m_renderProgram;
	GLuint m_VAO, m_VBO;
	GLint m_minCoordLoc;
	GLint m_maxCoordLoc;
	GLint m_pointSizeLoc;
	GLint m_pixelColorLoc;
};

