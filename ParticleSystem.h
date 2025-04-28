#pragma once

#include "Definition.h"

template <int DIM>
struct Particle {
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	MatrixX<DIM, DIM> affineMatrix;
	VectorX<DIM> position;
	VectorX<DIM> velocity;
	VectorX<DIM> buffer;
	Scalar radius;
	Scalar mass;
	bool isFixed;

	Particle()
		: affineMatrix(MatrixX<DIM, DIM>::Zero()), 
		position(VectorX<DIM>::Zero()), 
		velocity(VectorX<DIM>::Zero()),
		buffer(VectorX<DIM>::Zero()),
		radius(0.0), mass(1.0), isFixed(false) { }
};

template <int DIM>
class ParticleSystem {
public:
	ParticleSystem(Integer numParticles, const Scalar& radius = 0.01);

	void SpawnParticles(const VectorX<DIM>& minCoord, const VectorX<DIM>& maxCoord);

	inline const std::vector<Particle<DIM>>& GetParticles() const { return m_particle; }

	inline const VectorX<DIM>& GetPosition(Integer index) const { return m_particle[index].position; }
	inline const VectorX<DIM>& GetVelocity(Integer index) const { return m_particle[index].velocity; }
	inline const MatrixX<DIM, DIM>& GetAffineMatrix(Integer index) const { return m_particle[index].affineMatrix; }
	inline const Scalar& GetRadius(Integer index) const { return m_particle[index].radius; }
	inline const VectorX<DIM>& GetBuffer(Integer index) const { return m_particle[index].buffer; }
	inline const Particle<DIM>& GetParticle(Integer index) const { return m_particle[index]; }
	inline const Scalar& GetMass(Integer index) const { return m_particle[index].mass; }
	inline Particle<DIM>& GetParticle(Integer index) { return m_particle[index]; }

	inline void SetPosition(Integer index, const VectorX<DIM>& position) { m_particle[index].position = position; }
	inline void SetVelocity(Integer index, const VectorX<DIM>& velocity) { m_particle[index].velocity = velocity; }
	inline void SetAffineMatrix(Integer index, const MatrixX<DIM, DIM>& affineMatrix) { m_particle[index].affineMatrix = affineMatrix; }
	inline void SetRadius(Integer index, const Scalar& radius) { m_particle[index].radius = radius; }
	inline void SetBuffer(Integer index, const VectorX<DIM>& val) { m_particle[index].buffer = val; }
	inline void SetMass(Integer index, const Scalar& val) { m_particle[index].mass = val; }

	inline const Integer& NumParticles() const { return m_size; }

	void Resize(Integer size);

	void WriteToFile(const std::string& filename, bool isAppend = false) const;
	void LoadFromLine(const std::string& line);

private:
	Integer m_size;
	Scalar m_defaultRadius;
	std::vector<Particle<DIM>> m_particle;
};

