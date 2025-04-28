#pragma once

#include "Definition.h"
#include "Scene.h"
#include "ParticleSystem.h"
#include "HairState.h"
#include "Renderer.h"


template <int DIM>
class HairModel {
public:
	HairModel(const Scene<DIM>& scene);

	~HairModel();

	inline HairState<DIM>& GetHairState(Integer hairIndex) { return m_hairState[hairIndex]; }
	inline ParticleSystem<DIM>& GetParticleSystem() { return m_ps; }
	inline const HairParameter<DIM>& GetHairParameter() const { return m_params; }
	inline const Scene<DIM>& GetScene() const { return m_scene; }

	inline Particle<DIM>& GetParticle(Integer hairIndex, Integer vertexIndex) { return m_ps.GetParticle(m_hairOffset[hairIndex] + vertexIndex); }
	inline const VectorX<DIM>& GetPosition(Integer hairIndex, Integer vertexIndex) const { return m_ps.GetPosition(m_hairOffset[hairIndex] + vertexIndex); }
	inline const VectorX<DIM>& GetVelocity(Integer hairIndex, Integer vertexIndex) const { return m_ps.GetVelocity(m_hairOffset[hairIndex] + vertexIndex); }
	inline const Scalar& GetMass(Integer hairIndex, Integer vertexIndex) const { return m_ps.GetMass(m_hairOffset[hairIndex] + vertexIndex); }
	inline const Scalar& GetRadius(Integer hairIndex, Integer vertexIndex) const { return m_ps.GetRadius(m_hairOffset[hairIndex] + vertexIndex); }
	inline bool IsFixed(Integer hairIndex, Integer vertexIndex) const { return m_ps.GetParticle(m_hairOffset[hairIndex] + vertexIndex).isFixed; }
	inline void SetPosition(Integer hairIndex, Integer vertexIndex, const VectorX<DIM>& position) { m_ps.SetPosition(m_hairOffset[hairIndex] + vertexIndex, position); }
	inline void SetVelocity(Integer hairIndex, Integer vertexIndex, const VectorX<DIM>& velocity) { m_ps.SetVelocity(m_hairOffset[hairIndex] + vertexIndex, velocity); }
	inline void SetMass(Integer hairIndex, Integer vertexIndex, const Scalar& mass) { m_ps.SetMass(m_hairOffset[hairIndex] + vertexIndex, mass); }
	inline void SetRadius(Integer hairIndex, Integer vertexIndex, const Scalar& radius) { m_ps.SetRadius(m_hairOffset[hairIndex] + vertexIndex, radius); }

	inline Integer NumVertexOnHair(Integer hairIndex) const {
		return hairIndex >= m_params.numHairs - 1 ?
			m_ps.NumParticles() - m_hairOffset[hairIndex] :
			m_hairOffset[hairIndex + 1] - m_hairOffset[hairIndex];
	}

	inline Matrix2 GetRefRotationMatrix() const {
		Scalar c = std::cos(m_params.refRotation);
		Scalar s = std::sin(m_params.refRotation);

		Matrix2 mat;
		mat << c, -s,
			   s, c;
		return mat;
	}

	void Clear();

	void Init();

	void PreStep(const Scalar& dt);

	void Step(const Scalar& dt);

	void PostStep(const Scalar& dt);

	void Render() const;

private:
	const Scene<DIM>& m_scene;
	const HairParameter<DIM>& m_params;
	std::vector<HairState<DIM>> m_hairState;
	std::vector<Integer> m_hairOffset;	// offset of each hair in the particle system
	ParticleSystem<DIM> m_ps;

	// render data
	GLuint m_renderProgram;
	GLuint m_VAO, m_VBO, m_EBO;
	GLint m_minCoordLoc;
	GLint m_maxCoordLoc;
	GLint m_pointSizeLoc;
	GLint m_pixelColorLoc;
	Integer m_numRenderLines;
};

