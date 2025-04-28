#include "ParticleSystem.h"

template class ParticleSystem<2>;
template class ParticleSystem<3>;

template<int DIM>
ParticleSystem<DIM>::ParticleSystem(Integer numParticles, const Scalar& radius)
	: m_defaultRadius(2 * radius), m_size(numParticles){}

template<int DIM>
void ParticleSystem<DIM>::Resize(Integer size) {
	m_particle.resize(size);

	#pragma omp simd
	for (Integer index = m_size; index < size; ++index) {
		SetRadius(index, m_defaultRadius);
	}

	m_size = size;

}

template<>
void ParticleSystem<2>::SpawnParticles(const VectorX<2>& minCoord, const VectorX<2>& maxCoord) {
	Scalar interval = 0.5 * m_defaultRadius;
	Integer dx = (maxCoord.x() - minCoord.x()) / interval;
	Integer dy = (maxCoord.y() - minCoord.y()) / interval;
	Resize(min(dx * dy, m_size));

#ifdef _DEBUG
	std::cout << "Spawned " << m_size << " particles" << std::endl;
#endif

	for (Integer xi = 0; xi < dx; ++xi) {
		for (Integer yi = 0; yi < dy; ++yi) {
			Integer index = xi * dy + yi;
			if (index >= m_size) break;
			SetPosition(index, Vector2(xi, yi) * interval + minCoord);
			SetRadius(index, m_defaultRadius);
		}
	}
}

template<>
void ParticleSystem<3>::SpawnParticles(const VectorX<3>& minCoord, const VectorX<3>& maxCoord) {
	Scalar interval = 0.5 * m_defaultRadius;
	Vector3i dx = ((maxCoord - minCoord) / interval).cast<Integer>();
	Resize(min(m_size, dx.x() * dx.y() * dx.z()));

#ifdef _DEBUG
	std::cout << "Spawned " << m_size << " particles" << std::endl;
#endif

	for (Integer xi = 0; xi < dx.x(); ++xi) {
		for (Integer yi = 0; yi < dx.y(); ++yi) {
			for (Integer zi = 0; zi < dx.z(); ++zi) {
				Integer index = xi * dx.y() * dx.z() + yi * dx.z() + zi;
				if (index >= m_size) break;
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

	std::getline(ss, token, ',');
	for (size_t i = 0; i < m_size; ++i) {
		for (Integer di = 0; di < DIM; ++di) {
			if (!std::getline(ss, token, ',')) {
				continue;
				throw std::runtime_error("Failed to recognize line " + line);
			}
			m_particle[i].position[di] = std::stof(token);
		}
	}
}
