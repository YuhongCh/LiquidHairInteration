#pragma once

#include "Definition.h"
#include "ParticleSystem.h"
#include "Grid.h"
#include "MathUtils.h"


template <int DIM>
class HashNeighborSearch {
public:

	HashNeighborSearch(const ParticleSystem<DIM>& ps, const Grid<DIM>& grid, Integer scaleFactor = 2)
		: m_ps(ps), m_grid(grid), m_scaleFactor(scaleFactor), 
		m_size(scaleFactor* ps.NumParticles()), m_hashNeighbor(m_size), m_cellLock(m_size) {
		for (Integer index = 0; index < m_size; ++index) {
			omp_init_lock(&m_cellLock[index]);
		}
	}

	template <typename Func>
	void ForEachNeighborParticles(Integer pi, Func func) {
		VectorXi<DIM> minCoord, maxCoord;
		GetImpactRegion(pi, minCoord, maxCoord);

		VectorX<DIM> pos = m_ps.GetPosition(pi);
		Scalar r2 = m_ps.GetRadius(pi) * m_ps.GetRadius(pi);
		if constexpr (DIM == 2) {
			for (Integer xi = minCoord[0]; xi <= maxCoord[0]; ++xi) {
				for (Integer yi = minCoord[1]; yi <= maxCoord[1]; ++yi) {
					Integer slot = ((((p1 * xi) ^ (p2 * yi)) % m_size) + m_size) % m_size;
					for (Integer npi : m_hashNeighbor[slot]) {
						if (pi == npi) continue;
						if ((pos - m_ps.GetPosition(npi)).squaredNorm() >= r2) continue;
						func(pi, npi);
					}
				}
			}
		}
		else if constexpr (DIM == 3) {
			for (Integer xi = minCoord[0]; xi <= maxCoord[0]; ++xi) {
				for (Integer yi = minCoord[1]; yi <= maxCoord[1]; ++yi) {
					for (Integer zi = minCoord[2]; zi <= maxCoord[2]; ++zi) {
						Integer slot = ((((p1 * xi) ^ (p2 * yi) ^ (p3 * zi)) % m_size) + m_size) % m_size;
						for (Integer npi : m_hashNeighbor[slot]) {
							if (pi == npi) continue;
							if ((pos - m_ps.GetPosition(npi)).squaredNorm() >= r2) continue;
							func(pi, npi);
						}
					}
				}
			}
		}
		else {
			static_assert(DIM == 2 || DIM == 3, "Invalid template DIM value found, expected [2, 3].");
		}
	}

	void ReCompute() {
		if (m_size != m_scaleFactor * m_ps.NumParticles()) {
			for (auto& lk : m_cellLock) omp_destroy_lock(&lk);

			m_size = m_scaleFactor * m_ps.NumParticles();
			m_hashNeighbor.resize(m_size);
			m_cellLock.resize(m_size);

			for (Integer i = 0; i < m_size; ++i) {
				omp_init_lock(&m_cellLock[i]);
			}
		}

		for (Integer i = 0; i < m_size; ++i) m_hashNeighbor[i].clear();

		#pragma omp parallel for 
		for (Integer pi = 0; pi < m_ps.NumParticles(); ++pi) {
			VectorXi<DIM> minCoord, maxCoord;
			GetImpactRegion(pi, minCoord, maxCoord);

			if constexpr (DIM == 2) {
				for (Integer xi = minCoord[0]; xi <= maxCoord[0]; ++xi) {
					for (Integer yi = minCoord[1]; yi <= maxCoord[1]; ++yi) {
						Integer slot = ((((p1 * xi) ^ (p2 * yi)) % m_size) + m_size) % m_size;

						omp_set_lock(&m_cellLock[slot]);
						m_hashNeighbor[slot].push_back(pi);
						omp_unset_lock(&m_cellLock[slot]);
					}
				}
			}
			else if constexpr (DIM == 3) {
				for (Integer xi = minCoord[0]; xi <= maxCoord[0]; ++xi) {
					for (Integer yi = minCoord[1]; yi <= maxCoord[1]; ++yi) {
						for (Integer zi = minCoord[2]; zi <= maxCoord[2]; ++zi) {
							Integer slot = ((((p1 * xi) ^ (p2 * yi) ^ (p3 * zi)) % m_size) + m_size) % m_size;

							omp_set_lock(&m_cellLock[slot]);
							m_hashNeighbor[slot].push_back(pi);
							omp_unset_lock(&m_cellLock[slot]);
						}
					}
				}
			}
			else {
				static_assert(DIM == 2 || DIM == 3, "Invalid template DIM value found, expected [2, 3].");
			}

		}
	}

protected:
	void GetImpactRegion(Integer pi, VectorXi<DIM>& minCoord, VectorXi<DIM>& maxCoord) {
		VectorX<DIM> pos = m_ps.GetPosition(pi);
		VectorXi<DIM> coord = ((pos - m_grid.minCoord).cwiseQuotient(m_grid.dx)).cast<Integer>();

		for (Integer di = 0; di < DIM; ++di) {
			minCoord[di] = MathUtils::Clamp(coord[di] - 1, 0, m_grid.dimension[di] - 1);
			maxCoord[di] = MathUtils::Clamp(coord[di] + 1, 0, m_grid.dimension[di] - 1);
		}
		

		//Scalar radius = m_ps.GetRadius(pi);

		//for (Integer di = 0; di < DIM; ++di) {
		//	minCoord[di] = (Integer)std::floor((pos[di] - radius - m_grid.minCoord[di]) / m_grid.dx[di]);
		//	maxCoord[di] = (Integer)std::ceil((pos[di] + radius - m_grid.minCoord[di]) / m_grid.dx[di]);

		//	minCoord[di] = MathUtils::Clamp(minCoord[di], 0, m_grid.dimension[di] - 1);
		//	maxCoord[di] = MathUtils::Clamp(maxCoord[di], 0, m_grid.dimension[di] - 1);
		//}
	}

private:
	// use for hashing
	static inline constexpr Integer p1 = 73856093;
	static inline constexpr Integer p2 = 19349663;
	static inline constexpr Integer p3 = 83492791;

private:
	Integer m_scaleFactor;
	Integer m_size;
	const ParticleSystem<DIM>& m_ps;
	const Grid<DIM>& m_grid;
	std::vector<std::vector<Integer>> m_hashNeighbor;
	std::vector<omp_lock_t> m_cellLock;

};

