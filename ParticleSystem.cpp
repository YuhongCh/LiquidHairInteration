#include "ParticleSystem.h"

template class ParticleSystem<2>;
template class ParticleSystem<3>;

template<int DIM>
ParticleSystem<DIM>::ParticleSystem(int numParticles, const Scalar& radius)
	: m_size(numParticles), m_defaultRadius(2 * radius){

	m_particle.resize(numParticles);
}

template<int DIM>
void ParticleSystem<DIM>::Resize(Integer size) {
	m_particle.resize(size);
	m_size = size;
}

template<>
void ParticleSystem<2>::SpawnParticles(const VectorX<2>& minCoord, const VectorX<2>& maxCoord) {
	Scalar interval = 0.5 * m_defaultRadius;
	Integer dx = (maxCoord.x() - minCoord.x()) / interval;
	Integer dy = (maxCoord.y() - minCoord.y()) / interval;
	Resize(dx * dy);

#ifdef _DEBUG
	std::cout << "Spawned " << m_size << " particles" << std::endl;
#endif

	for (Integer xi = 0; xi < dx; ++xi) {
		for (Integer yi = 0; yi < dy; ++yi) {
			Integer index = xi * dy + yi;
			SetPosition(index, Vector2(xi, yi) * interval + minCoord);
			SetRadius(index, m_defaultRadius);
		}
	}
}

template<>
void ParticleSystem<3>::SpawnParticles(const VectorX<3>& minCoord, const VectorX<3>& maxCoord) {
	Scalar interval = 0.5 * m_defaultRadius;
	Vector3i dx = ((maxCoord - minCoord) / interval).cast<Integer>();
	Resize(dx.x() * dx.y() * dx.z());

#ifdef _DEBUG
	std::cout << "Spawned " << m_size << " particles" << std::endl;
#endif

	for (Integer xi = 0; xi < dx.x(); ++xi) {
		for (Integer yi = 0; yi < dx.y(); ++yi) {
			for (Integer zi = 0; zi < dx.z(); ++zi) {
				Integer index = xi * dx.y() * dx.z() + yi * dx.z() + zi;
				SetPosition(index, Vector3(xi, yi, zi) * interval + minCoord);
				SetRadius(index, m_defaultRadius);
			}
		}
	}
}

template<int DIM>
void ParticleSystem<DIM>::WriteToFile(const std::string& filename, bool isAppend) const {
	std::ofstream writer;
	if (isAppend) writer.open(filename, std::ios::out | std::ios::app);
	else writer.open(filename, std::ios::out | std::ios::trunc);

	if (!writer.is_open()) {
		throw std::runtime_error("Failed to open file: " + filename);
	}

	writer << m_particle.size() << "\n";
	for (const Particle<DIM>& particle : m_particle) {
		for (Integer di = 0; di < DIM; ++di) {
			writer << particle.position[di] << ",";
		}
	}
	writer << "\n";
	writer.close();
}

template<int DIM>
void ParticleSystem<DIM>::LoadFromLine(const std::string& line) {
	std::string token;
	std::istringstream ss(line);

	size_t numParticles = std::stoul(token);
	if (numParticles != m_size) Resize(numParticles);

	std::getline(ss, token, ',');
	for (size_t i = 0; i < numParticles; ++i) {
		for (Integer di = 0; di < DIM; ++di) {
			if (!std::getline(ss, token, ',')) {
				throw std::runtime_error("Failed to recognize line " + line);
			}
			m_particle[i].position[di] = std::stof(token);
		}
	}
}
