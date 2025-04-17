#pragma once

#include "Definition.h"
#include "ParticleSystem.h"
#include "Grid.h"
#include "MathUtils.h"


template <int DIM>
class HashNeighborSearch {
public:

	HashNeighborSearch(const ParticleSystem<DIM>& ps, const Grid<DIM>& grid, Integer scaleFactor = 3)
		: m_ps(ps), m_grid(grid), m_scaleFactor(scaleFactor), 
		m_size(scaleFactor* ps.NumParticles()), m_hashNeighbor(scaleFactor * ps.NumParticles()){ }

	Integer Hash(const VectorXi<DIM>& coord) const;

	template <typename Func>
	void ForEachNeighborParticles(Integer pi, Func func) {
		VectorXi<DIM> minCoord, maxCoord;
		GetImpactRegion(pi, minCoord, maxCoord);

		if constexpr (DIM == 2) {
			for (Integer xi = minCoord[0]; xi <= maxCoord[0]; ++xi) {
				for (Integer yi = minCoord[1]; yi <= maxCoord[1]; ++yi) {
					Integer slot = Hash(Vector2i(xi, yi));
					for (Integer npi : m_hashNeighbor[slot]) {
						if (pi == npi) continue;
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

	void MoveParticle(Integer pi, const VectorX<DIM>& from, const VectorX<DIM>& to) {
		VectorXi<DIM> minCoordFrom, maxCoordFrom;
		VectorXi<DIM> minCoordTo, maxCoordTo;

		Scalar radius = m_ps.GetRadius(pi);
		for (Integer di = 0; di < DIM; ++di) {
			minCoordFrom[di] = (Integer)std::floor((from[di] - radius - m_grid.minCoord[di]) / m_grid.dx[di]);
			maxCoordFrom[di] = (Integer)std::ceil((from[di] + radius - m_grid.minCoord[di]) / m_grid.dx[di]);
			minCoordFrom[di] = MathUtils::Clamp(minCoordFrom[di], 0, m_grid.dimension[di] - 1);
			maxCoordFrom[di] = MathUtils::Clamp(maxCoordFrom[di], 0, m_grid.dimension[di] - 1);

			minCoordTo[di] = (Integer)std::floor((to[di] - radius - m_grid.minCoord[di]) / m_grid.dx[di]);
			maxCoordTo[di] = (Integer)std::ceil((to[di] + radius - m_grid.minCoord[di]) / m_grid.dx[di]);
			minCoordTo[di] = MathUtils::Clamp(minCoordTo[di], 0, m_grid.dimension[di] - 1);
			maxCoordTo[di] = MathUtils::Clamp(maxCoordTo[di], 0, m_grid.dimension[di] - 1);
		}

		if constexpr (DIM == 2) {
			for (Integer xi = minCoordFrom[0]; xi <= maxCoordFrom[0]; ++xi) {
				for (Integer yi = minCoordFrom[1]; yi <= maxCoordFrom[1]; ++yi) {
					Integer slot = Hash(Vector2i(xi, yi));
					for (Integer i = 0; i < m_hashNeighbor[slot].size(); ++i) {
						if (m_hashNeighbor[slot][i] == pi) {
							std::swap(m_hashNeighbor[slot][i], m_hashNeighbor[slot].back());
							m_hashNeighbor[slot].pop_back();
							break;
						}
					}
				}
			}

			for (Integer xi = minCoordTo[0]; xi <= maxCoordTo[0]; ++xi) {
				for (Integer yi = minCoordTo[1]; yi <= maxCoordTo[1]; ++yi) {
					Integer slot = Hash(Vector2i(xi, yi));
					m_hashNeighbor[slot].push_back(pi);
				}
			}
		}
		else if constexpr (DIM == 3) {
			for (Integer xi = minCoordFrom[0]; xi <= maxCoordFrom[0]; ++xi) {
				for (Integer yi = minCoordFrom[1]; yi <= maxCoordFrom[1]; ++yi) {
					for (Integer zi = minCoordFrom[2]; zi <= maxCoordFrom[2]; ++zi) {
						Integer slot = Hash(Vector3i(xi, yi, zi));
						for (Integer i = 0; i < m_hashNeighbor[slot].size(); ++i) {
							if (m_hashNeighbor[slot][i] == pi) {
								std::swap(m_hashNeighbor[slot][i], m_hashNeighbor[slot].back());
								m_hashNeighbor[slot].pop_back();
								break;
							}
						}
					}
				}
			}

			for (Integer xi = minCoordTo[0]; xi <= maxCoordTo[0]; ++xi) {
				for (Integer yi = minCoordTo[1]; yi <= maxCoordTo[1]; ++yi) {
					for (Integer zi = minCoordTo[2]; zi <= maxCoordTo[2]; ++zi) {
						Integer slot = Hash(Vector3i(xi, yi, zi));
						m_hashNeighbor[slot].push_back(pi);
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
			m_size = m_scaleFactor * m_ps.NumParticles();
			m_hashNeighbor.resize(m_size);
		}

		for (Integer i = 0; i < m_size; ++i) m_hashNeighbor[i].clear();

		#pragma omp parallel for 
		for (Integer pi = 0; pi < m_ps.NumParticles(); ++pi) {
			VectorXi<DIM> minCoord, maxCoord;
			GetImpactRegion(pi, minCoord, maxCoord);

			if constexpr (DIM == 2) {
				for (Integer xi = minCoord[0]; xi <= maxCoord[0]; ++xi) {
					for (Integer yi = minCoord[1]; yi <= maxCoord[1]; ++yi) {
						Integer slot = Hash(Vector2i(xi, yi));

						#pragma omp critical
						m_hashNeighbor[slot].push_back(pi);
					}
				}
			}
			else if constexpr (DIM == 3) {
				for (Integer xi = minCoord[0]; xi <= maxCoord[0]; ++xi) {
					for (Integer yi = minCoord[1]; yi <= maxCoord[1]; ++yi) {
						for (Integer zi = minCoord[2]; zi <= maxCoord[2]; ++zi) {
							Integer slot = ((((p1 * xi) ^ (p2 * yi) ^ (p3 * zi)) % m_size) + m_size) % m_size;

							#pragma omp critical
							m_hashNeighbor[slot].push_back(pi);
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
		Scalar radius = m_ps.GetRadius(pi);

		for (Integer di = 0; di < DIM; ++di) {
			minCoord[di] = (Integer)std::floor((pos[di] - radius - m_grid.minCoord[di]) / m_grid.dx[di]);
			maxCoord[di] = (Integer)std::ceil((pos[di] + radius - m_grid.minCoord[di]) / m_grid.dx[di]);

			minCoord[di] = MathUtils::Clamp(minCoord[di], 0, m_grid.dimension[di] - 1);
			maxCoord[di] = MathUtils::Clamp(maxCoord[di], 0, m_grid.dimension[di] - 1);
		}
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

};

