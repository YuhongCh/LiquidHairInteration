#include "PressureSolver.h"
#include "MathUtils.h"


#pragma region Jacobi Solver
template class JacobiSolver<2>;

#pragma region Jacobi Solver 2D

void JacobiSolver<2>::BuildRHS() {
	for (Integer xi = 0; xi < m_grid.dimension.x(); ++xi) {
		for (Integer yi = 0; yi < m_grid.dimension.y(); ++yi) {
			m_rhs[xi][yi] = 0.0;
			if (m_grid.liquidPhi[xi][yi] > 0.0) continue;

			m_rhs[xi][yi] -= m_grid.liquidWx[xi + 1][yi] * m_grid.vx[xi + 1][yi] / m_grid.dx.x();
			m_rhs[xi][yi] += m_grid.liquidWx[xi][yi] * m_grid.vx[xi][yi] / m_grid.dx.x();
			m_rhs[xi][yi] -= m_grid.liquidWy[xi][yi + 1] * m_grid.vy[xi][yi + 1] / m_grid.dx.y();
			m_rhs[xi][yi] += m_grid.liquidWy[xi][yi] * m_grid.vy[xi][yi] / m_grid.dx.y();
		}
	}
}

void JacobiSolver<2>::Solve(const Scalar& liquidDensity, const Scalar& dt, bool warmStart) {
	BuildRHS();

	if (!warmStart) {
		for (Integer xi = 0; xi < m_grid.dimension.x(); ++xi) {
			m_buffer0[xi].assign(m_grid.dimension.y(), 0.0);
		}
	}

	for (Integer iter = 0; iter < m_maxIteration; ++iter) {
		#pragma omp parallel for
		for (Integer xi = 1; xi < m_grid.dimension.x() - 1; ++xi) {
			for (Integer yi = 1; yi < m_grid.dimension.y() - 1; ++yi) {
				Scalar leftValue = 0.0, rightValue = 0.0, topValue = 0.0, botValue = 0.0, diagValue = 0.0;
				Scalar centerPhi = m_grid.liquidPhi[xi][yi];
				if (centerPhi > 0.0) continue;

				// right
				if (m_grid.liquidWx[xi + 1][yi] > 0.0) {
					Scalar term = m_grid.liquidWx[xi + 1][yi] * dt / (liquidDensity * m_grid.dx.x() * m_grid.dx.x());
					Scalar phi = m_grid.liquidPhi[xi + 1][yi];
					if (phi < 0.0) {
						rightValue -= term * m_buffer0[xi + 1][yi];
						diagValue += term;
					}
					else {
						Scalar frac = MathUtils::FractionInRegion(centerPhi, phi);
						diagValue += term / frac;
					}
				}

				// left
				if (m_grid.liquidWx[xi][yi] > 0.0) {
					Scalar term = m_grid.liquidWx[xi][yi] * dt / (liquidDensity * m_grid.dx.x() * m_grid.dx.x());
					Scalar phi = m_grid.liquidPhi[xi - 1][yi];
					if (phi < 0.0) {
						leftValue -= term * m_buffer0[xi - 1][yi];
						diagValue += term;
					}
					else {
						Scalar frac = MathUtils::FractionInRegion(centerPhi, phi);
						diagValue += term / frac;
					}
				}

				// top
				if (m_grid.liquidWy[xi][yi + 1] > 0.0) {
					Scalar term = m_grid.liquidWy[xi][yi + 1] * dt / (liquidDensity * m_grid.dx.y() * m_grid.dx.y());
					Scalar phi = m_grid.liquidPhi[xi][yi + 1];
					if (phi < 0.0) {
						topValue -= term * m_buffer0[xi][yi + 1];
						diagValue += term;
					}
					else {
						Scalar frac = MathUtils::FractionInRegion(centerPhi, phi);
						diagValue += term / frac;
					}
				}

				// bottom
				if (m_grid.liquidWy[xi][yi] > 0.0) {
					Scalar term = m_grid.liquidWy[xi][yi] * dt / (liquidDensity * m_grid.dx.y() * m_grid.dx.y());
					Scalar phi = m_grid.liquidPhi[xi][yi - 1];
					if (phi < 0.0) {
						botValue -= term * m_buffer0[xi][yi - 1];
						diagValue += term;
					}
					else {
						Scalar frac = MathUtils::FractionInRegion(centerPhi, phi);
						diagValue += term / frac;
					}
				}

				if (!MathUtils::IsSmall(std::abs(diagValue))) {
					m_buffer1[xi][yi] = (1.0 - m_relaxFactor) * m_buffer0[xi][yi] + m_relaxFactor * (m_rhs[xi][yi] - (leftValue + rightValue + topValue + botValue)) / diagValue;
				}
			}
		}
		std::swap(m_buffer0, m_buffer1);
	}

	for (Integer xi = 0; xi < m_grid.dimension.x(); ++xi) {
		for (Integer yi = 0; yi < m_grid.dimension.y(); ++yi) {
			m_grid.pressure[xi][yi] = m_buffer0[xi][yi];
		}
	}
}

#pragma endregion
#pragma region Jacobi Solver 3D

#pragma endregion

#pragma endregion


#pragma region Conjugate Gradient Solver
template class CGSolver<2>;
template class CGSolver<3>;

#pragma region Conjugate Gradient Solver 2D

CGSolver<2>::CGSolver(Grid<2>& grid, Integer maxIteration,
					  Scalar tolerance, PrecondType type,
					  Integer maxLevel, Integer numPreSmooth, Integer numPostSmooth, Integer numFinalSmooth)
	: m_grid(grid), m_maxIteration(maxIteration), m_tolerance(tolerance), m_type(type),
	  m_maxLevel(max(1, maxLevel)), m_numPreSmooth(numPreSmooth), m_numPostSmooth(numPostSmooth), m_numFinalSmooth(numFinalSmooth) {

	m_rhs = std::vector<std::vector<Scalar>>(m_grid.dimension.x(), std::vector<Scalar>(m_grid.dimension.y()));
	m_p = std::vector<std::vector<Scalar>>(m_grid.dimension.x(), std::vector<Scalar>(m_grid.dimension.y()));
	m_x = std::vector<std::vector<Scalar>>(m_grid.dimension.x(), std::vector<Scalar>(m_grid.dimension.y()));

	m_maxLevel = (m_type != PrecondType::Multigrid) ? 1 : m_maxLevel;
	m_r.resize(m_maxLevel);
	m_z.resize(m_maxLevel);
	m_Adiag.resize(m_maxLevel);
	m_Acoef.resize(m_maxLevel);

	for (Integer level = 0; level < m_maxLevel; ++level) {
		Integer width = m_grid.dimension.x() / (1 << level);
		Integer height = m_grid.dimension.y() / (1 << level);
		m_r[level] = std::vector<std::vector<Scalar>>(width, std::vector<Scalar>(height));
		m_z[level] = std::vector<std::vector<Scalar>>(width, std::vector<Scalar>(height));
		m_Adiag[level] = std::vector<std::vector<Scalar>>(width, std::vector<Scalar>(height));
		m_Acoef[level] = std::vector<std::vector<Vector4>>(width, std::vector<Vector4>(height));
	}

	m_oldRZ = 0.0;
	m_newRZ = 0.0;
	m_alpha = 0.0;
	m_beta = 0.0;
}

Scalar CGSolver<2>::Reduce(const std::vector<std::vector<Scalar>>& vec0, const std::vector<std::vector<Scalar>>& vec1) const {
	Scalar sum = 0.0;
	#pragma omp parallel for reduction(+:sum)
	for (Integer xi = 0; xi < vec0.size(); ++xi) {
		for (Integer yi = 0; yi < vec0[xi].size(); ++yi) {
			sum += vec0[xi][yi] * vec1[xi][yi];
		}
	}
	return sum;
}

Scalar CGSolver<2>::GetApproxLiquidPhi(Integer xi, Integer yi, Integer level) const {
	if (level < 0) throw std::invalid_argument("Negative level is not accepted");
	if (level == 0) return m_grid.liquidPhi[xi][yi];
	Integer factor = 1 << level;
	Integer baseX = factor * xi + (1 << (level - 1)) - 1, baseY = factor * yi + (1 << (level - 1)) - 1;
	return MathUtils::Bilerp(m_grid.liquidPhi[baseX][baseY], m_grid.liquidPhi[baseX + 1][baseY],
							 m_grid.liquidPhi[baseX][baseY + 1], m_grid.liquidPhi[baseX + 1][baseY + 1], 0.5, 0.5);
}

Scalar CGSolver<2>::GetApproxWx(Integer xi, Integer yi, Integer level) const {
	if (level < 0) throw std::invalid_argument("Negative level is not accepted");
	if (level == 0) return m_grid.liquidWx[xi][yi];
	Integer factor = 1 << level;
	Integer baseX = factor * xi, baseY = factor * yi + (1 << (level - 1)) - 1;
	return 0.5 * (m_grid.liquidWx[baseX][baseY] + m_grid.liquidWx[baseX][baseY + 1]);
}

Scalar CGSolver<2>::GetApproxWy(Integer xi, Integer yi, Integer level) const {
	if (level < 0) throw std::invalid_argument("Negative level is not accepted");
	if (level == 0) return m_grid.liquidWy[xi][yi];
	Integer factor = 1 << level;
	Integer baseX = factor * xi + (1 << (level - 1)) - 1, baseY = factor * yi;
	return 0.5 * (m_grid.liquidWy[baseX][baseY] + m_grid.liquidWy[baseX + 1][baseY]);
}

void CGSolver<2>::BuildLHS(const Scalar& liquidDensity, const Scalar& dt) {
	// build levels
	for (Integer level = 0; level < m_maxLevel; ++level) {
		Integer factor = 1 << level;
		Vector2i dimension(m_grid.dimension.x() / factor, m_grid.dimension.y() / factor);
		Vector2 dx = factor * m_grid.dx;

		#pragma omp parallel for
		for (Integer xi = 0; xi < dimension.x(); ++xi) {
			for (Integer yi = 0; yi < dimension.y(); ++yi) {
				Scalar centerPhi = GetApproxLiquidPhi(xi, yi, level);
				if (centerPhi > 0.0) continue;

				// right
				if (xi + 1 < dimension.x()) {
					Scalar rightWx = GetApproxWx(xi + 1, yi, level);
					if (rightWx > 0.0) {
						Scalar term = rightWx * dt / (liquidDensity * dx.x() * dx.x());
						Scalar phi = GetApproxLiquidPhi(xi + 1, yi, level);
						if (phi < 0.0) {
							m_Acoef[level][xi][yi][0] -= term;
							m_Adiag[level][xi][yi] += term;
						}
						else {
							Scalar frac = max(0.1, MathUtils::FractionInRegion(centerPhi, phi));
							m_Adiag[level][xi][yi] += term / frac;
						}
					}	
				}

				// left
				if (xi - 1 >= 0) {
					Scalar leftWx = GetApproxWx(xi, yi, level);
					if (leftWx > 0.0) {
						Scalar term = leftWx * dt / (liquidDensity * dx.x() * dx.x());
						Scalar phi = GetApproxLiquidPhi(xi - 1, yi, level);
						if (phi < 0.0) {
							m_Acoef[level][xi][yi][1] -= term;
							m_Adiag[level][xi][yi] += term;
						}
						else {
							Scalar frac = max(0.1, MathUtils::FractionInRegion(centerPhi, phi));
							m_Adiag[level][xi][yi] += term / frac;
						}
					}
				}

				// top
				if (yi + 1 < dimension.y()) {
					Scalar topWy = GetApproxWy(xi, yi + 1, level);
					if (topWy > 0.0) {
						Scalar term = topWy * dt / (liquidDensity * dx.y() * dx.y());
						Scalar phi = GetApproxLiquidPhi(xi, yi + 1, level);
						if (phi < 0.0) {
							m_Acoef[level][xi][yi][2] -= term;
							m_Adiag[level][xi][yi] += term;
						}
						else {
							Scalar frac = max(0.1, MathUtils::FractionInRegion(centerPhi, phi));
							m_Adiag[level][xi][yi] += term / frac;
						}
					}
				}

				// bottom
				if (yi - 1 >= 0) {
					Scalar botWy = GetApproxWy(xi, yi, level);
					if (botWy > 0.0) {
						Scalar term = botWy * dt / (liquidDensity * dx.y() * dx.y());
						Scalar phi = GetApproxLiquidPhi(xi, yi - 1, level);
						if (phi < 0.0) {
							m_Acoef[level][xi][yi][3] -= term;
							m_Adiag[level][xi][yi] += term;
						}
						else {
							Scalar frac = max(0.1, MathUtils::FractionInRegion(centerPhi, phi));
							m_Adiag[level][xi][yi] += term / frac;
						}
					}	
				}
			}
		}

	}
}

void CGSolver<2>::Init() {
	m_oldRZ = m_newRZ = m_alpha = m_beta = 0.0;

	for (Integer xi = 0; xi < m_grid.dimension.x(); ++xi) {
		m_p[xi].assign(m_grid.dimension.y(), Scalar(0.0));
		m_x[xi].assign(m_grid.dimension.y(), Scalar(0.0));
		m_rhs[xi].assign(m_grid.dimension.y(), Scalar(0.0));
	}

	#pragma omp parallel for
	for (Integer level = 0; level < m_maxLevel; ++level) {
		Integer width = m_grid.dimension.x() / (1 << level);
		Integer height = m_grid.dimension.y() / (1 << level);

		for (Integer xi = 0; xi < width; ++xi) {
			m_r[level][xi].assign(height, 0.0);
			m_z[level][xi].assign(height, 0.0);
			m_Adiag[level][xi].assign(height, 0.0);
			m_Acoef[level][xi].assign(height, Vector4::Zero());
		}
	}
}

void CGSolver<2>::BuildRHS() {
	#pragma omp parallel for
	for (Integer xi = 0; xi < m_grid.dimension.x(); ++xi) {
		for (Integer yi = 0; yi < m_grid.dimension.y(); ++yi) {
			m_rhs[xi][yi] = 0.0;
			if (m_grid.liquidPhi[xi][yi] > 0.0) continue;
			m_rhs[xi][yi] -= m_grid.liquidWx[xi + 1][yi] * m_grid.vx[xi + 1][yi] / m_grid.dx.x();
			m_rhs[xi][yi] += m_grid.liquidWx[xi][yi] * m_grid.vx[xi][yi] / m_grid.dx.x();
			m_rhs[xi][yi] -= m_grid.liquidWy[xi][yi + 1] * m_grid.vy[xi][yi + 1] / m_grid.dx.y();
			m_rhs[xi][yi] += m_grid.liquidWy[xi][yi] * m_grid.vy[xi][yi] / m_grid.dx.y();
		}
	}
}

Scalar CGSolver<2>::GetAx(const std::vector<std::vector<Scalar>>& grid, Integer xi, Integer yi, Integer level) const {
	Integer width = m_grid.dimension.x() / (1 << level);
	Integer height = m_grid.dimension.y() / (1 << level);
#ifdef _DEBUG
	if (xi < 0 || yi < 0 || xi >= width || yi >= height || level < 0 || level >= m_maxLevel) {
		throw std::invalid_argument("Get invalid or unexpected argument here");
	}

#endif
	const Vector4& coef = m_Acoef[level][xi][yi];
	Scalar result = grid[xi][yi] * m_Adiag[level][xi][yi];
	if (xi + 1 < width) result += grid[xi + 1][yi] * coef[0];
	if (xi - 1 >= 0) result += grid[xi - 1][yi] * coef[1];
	if (yi + 1 < height) result += grid[xi][yi + 1] * coef[2];
	if (yi - 1 >= 0) result += grid[xi][yi - 1] * coef[3];
	return result;
}

void CGSolver<2>::ApplyPreconditioner() {
	switch (m_type) {
		case PrecondType::None: {
			for (Integer xi = 0; xi < m_grid.dimension.x(); ++xi) {
				for (Integer yi = 0; yi < m_grid.dimension.y(); ++yi) {
					if (m_grid.liquidPhi[xi][yi] > 0.0) continue;
					m_z[0][xi][yi] = m_r[0][xi][yi];
				}
			}
			break;
		}
		case PrecondType::Jacobi: {
			for (Integer xi = 0; xi < m_grid.dimension.x(); ++xi) {
				for (Integer yi = 0; yi < m_grid.dimension.y(); ++yi) {
					Scalar denom = m_Adiag[0][xi][yi];
					if (MathUtils::IsSmall(denom)) continue;
					m_z[0][xi][yi] = m_r[0][xi][yi] / denom;
				}
			}
			break;
		}
		case PrecondType::Multigrid: {
			for (Integer level = 1; level < m_maxLevel; ++level) {
				for (Integer iter = 0; iter < m_numPreSmooth; ++iter) {
					Smooth(level - 1, 0);
					Smooth(level - 1, 1);
				}

				Restrict(level);
			}

			for (Integer iter = 0; iter < m_numFinalSmooth; ++iter) {
				Smooth(m_maxLevel - 1, 0);
				Smooth(m_maxLevel - 1, 1);
			}

			for (Integer level = m_maxLevel - 2; level >= 0; --level) {
				Prolong(level);
				for (Integer iter = 0; iter < m_numPostSmooth; ++iter) {
					Smooth(level, 0);
					Smooth(level, 1);
				}
			}
			break;
		}
		default: {
			throw std::runtime_error("Failed to recognize the preconditioner type");
		}
	}
}

void  CGSolver<2>::Prolong(Integer toLevel) {
	Integer fromLevel = toLevel + 1;
	Integer width = m_grid.dimension.x() / (1 << fromLevel);
	Integer height = m_grid.dimension.y() / (1 << fromLevel);

	#pragma omp parallel for 
	for (Integer xi = 0; xi < width; ++xi) {
		for (Integer yi = 0; yi < height; ++yi) {
			m_z[toLevel][2 * xi][2 * yi] += m_z[fromLevel][xi][yi];
			m_z[toLevel][2 * xi + 1][2 * yi] += m_z[fromLevel][xi][yi];
			m_z[toLevel][2 * xi][2 * yi + 1] += m_z[fromLevel][xi][yi];
			m_z[toLevel][2 * xi + 1][2 * yi + 1] += m_z[fromLevel][xi][yi];
		}
	}
}

void CGSolver<2>::Restrict(Integer toLevel) {
	Integer fromLevel = toLevel - 1;
	Integer width = m_grid.dimension.x() / (1 << toLevel);
	Integer height = m_grid.dimension.y() / (1 << toLevel);
	#pragma omp parallel for 
	for (Integer xi = 0; xi < width; ++xi) {
		for (Integer yi = 0; yi < height; ++yi) {
			m_r[toLevel][xi][yi] = 0.0;
			m_z[toLevel][xi][yi] = 0.0;

			m_r[toLevel][xi][yi] += 0.25 * (m_r[fromLevel][2 * xi][2 * yi] - GetAx(m_z[fromLevel], 2 * xi, 2 * yi, fromLevel));
			m_r[toLevel][xi][yi] += 0.25 * (m_r[fromLevel][2 * xi + 1][2 * yi] - GetAx(m_z[fromLevel], 2 * xi + 1, 2 * yi, fromLevel));
			m_r[toLevel][xi][yi] += 0.25 * (m_r[fromLevel][2 * xi][2 * yi + 1] - GetAx(m_z[fromLevel], 2 * xi, 2 * yi + 1, fromLevel));
			m_r[toLevel][xi][yi] += 0.25 * (m_r[fromLevel][2 * xi + 1][2 * yi + 1] - GetAx(m_z[fromLevel], 2 * xi + 1, 2 * yi + 1, fromLevel));
		}
	}
}

void CGSolver<2>::Smooth(Integer level, Integer phase) {
	Integer width = m_grid.dimension.x() / (1 << level);
	Integer height = m_grid.dimension.y() / (1 << level);

	#pragma omp parallel for
	for (Integer xi = 0; xi < width; ++xi) {
		for (Integer yi = 0; yi < height; ++yi) {
			if (((xi + yi) & 0b1) == phase) {
				if (!MathUtils::IsSmall(m_Adiag[level][xi][yi])) {
					Scalar term = (m_r[level][xi][yi] - GetAx(m_z[level], xi, yi, level)) / m_Adiag[level][xi][yi];
					m_z[level][xi][yi] = 0.5 * m_z[level][xi][yi] + 0.5 * term;
				}
			}
		}
	}
}

void CGSolver<2>::Solve(const Scalar& liquidDensity, const Scalar& dt) {
	// Initialize
	Init();
	BuildRHS();
	BuildLHS(liquidDensity, dt);

	Integer numLiquidCell = 0;
	for (Integer xi = 0; xi < m_grid.dimension.x(); ++xi) {
		for (Integer yi = 0; yi < m_grid.dimension.y(); ++yi) {
			if (m_grid.liquidPhi[xi][yi] > 0.0) continue;
			m_r[0][xi][yi] = m_rhs[xi][yi];
			++numLiquidCell;
		}
	}
	ApplyPreconditioner();
	for (Integer xi = 0; xi < m_grid.dimension.x(); ++xi) {
		for (Integer yi = 0; yi < m_grid.dimension.y(); ++yi) {
			if (m_grid.liquidPhi[xi][yi] > 0.0) continue;
			m_p[xi][yi] = m_z[0][xi][yi];
		}
	}

	m_oldRZ = Reduce(m_r[0], m_z[0]);
	Scalar residual =  std::sqrt(Reduce(m_r[0], m_r[0])) / numLiquidCell;

	// start real solving
	for (Integer iter = 0; iter < m_maxIteration; ++iter) {
		if (residual < m_tolerance) break;

		Scalar denom = 0.0;
		for (Integer xi = 1; xi < m_grid.dimension.x() - 1; ++xi) {
			for (Integer yi = 1; yi < m_grid.dimension.y() - 1; ++yi) {
				if (m_grid.liquidPhi[xi][yi] > 0.0) continue;
				denom += m_p[xi][yi] * GetAx(m_p, xi, yi, 0);
			}
		}
		m_alpha = m_oldRZ / denom;
#ifdef _DEBUG
		ASSERT_MSG(!std::isnan(m_alpha) && !std::isinf(m_alpha), 
			"Receive alpha=" + std::to_string(m_alpha) + " with oldRZ=" + std::to_string(m_oldRZ) + " and denom=" + std::to_string(denom));
#endif

		numLiquidCell = 0;
		residual = 0.0;
		#pragma omp parallel for  schedule(dynamic)
		for (Integer xi = 1; xi < m_grid.dimension.x() - 1; ++xi) {
			for (Integer yi = 1; yi < m_grid.dimension.y() - 1; ++yi) {
				if (m_grid.liquidPhi[xi][yi] > 0.0) continue;

				m_x[xi][yi] += m_alpha * m_p[xi][yi];
				m_r[0][xi][yi] -= m_alpha * GetAx(m_p, xi, yi, 0);
				residual += m_r[0][xi][yi] * m_r[0][xi][yi];
				++numLiquidCell;
			}
		}
		residual = std::sqrt(residual) / numLiquidCell;

		ApplyPreconditioner();
		m_newRZ = Reduce(m_r[0], m_z[0]);
		m_beta = m_newRZ / m_oldRZ;

#ifdef _DEBUG
		ASSERT_MSG(!std::isnan(m_beta) && !std::isinf(m_beta),
			"Receive beta=" + std::to_string(m_beta) + " with oldRZ=" + std::to_string(m_oldRZ) + " and newRZ=" + std::to_string(m_newRZ));
#endif

		#pragma omp parallel for  schedule(dynamic)
		for (Integer xi = 0; xi < m_grid.dimension.x(); ++xi) {
			for (Integer yi = 0; yi < m_grid.dimension.y(); ++yi) {
				if (m_grid.liquidPhi[xi][yi] > 0.0) continue;
				m_p[xi][yi] = m_z[0][xi][yi] + m_beta * m_p[xi][yi];
			}
		}
		m_oldRZ = m_newRZ;

#ifdef _DEBUG
		std::cout << "Iteration " << iter << ": residual=" << residual << std::endl;
#endif
	}

#ifdef _DEBUG
	ASSERT_MSG(residual < m_tolerance, "Get too large the residual: " + std::to_string(residual));
#endif

	for (Integer xi = 0; xi < m_grid.dimension.x(); ++xi) {
		for (Integer yi = 0; yi < m_grid.dimension.y(); ++yi) {
			m_grid.pressure[xi][yi] = m_x[xi][yi];
		}
	}
}

#pragma endregion

#pragma region Conjugate Gradient Solver 3D

CGSolver<3>::CGSolver(Grid<3>& grid, Integer maxIteration,
	Scalar tolerance, PrecondType type,
	Integer maxLevel, Integer numPreSmooth, Integer numPostSmooth, Integer numFinalSmooth)
	: m_grid(grid), m_maxIteration(maxIteration), m_tolerance(tolerance), m_type(type),
	m_maxLevel(max(1, maxLevel)), m_numPreSmooth(numPreSmooth), m_numPostSmooth(numPostSmooth), m_numFinalSmooth(numFinalSmooth) {

	m_rhs = std::vector<std::vector<std::vector<Scalar>>>(m_grid.dimension.x(), std::vector<std::vector<Scalar>>(m_grid.dimension.y(), std::vector<Scalar>(m_grid.dimension.z())));
	m_p = std::vector<std::vector<std::vector<Scalar>>>(m_grid.dimension.x(), std::vector<std::vector<Scalar>>(m_grid.dimension.y(), std::vector<Scalar>(m_grid.dimension.z())));
	m_x = std::vector<std::vector<std::vector<Scalar>>>(m_grid.dimension.x(), std::vector<std::vector<Scalar>>(m_grid.dimension.y(), std::vector<Scalar>(m_grid.dimension.z())));

	m_maxLevel = (m_type != PrecondType::Multigrid) ? 1 : m_maxLevel;
	m_r.resize(m_maxLevel);
	m_z.resize(m_maxLevel);
	m_Adiag.resize(m_maxLevel);
	m_Acoef.resize(m_maxLevel);

	for (Integer level = 0; level < m_maxLevel; ++level) {
		Vector3i dimension = m_grid.dimension / (1 << level);
		m_r[level] = std::vector<std::vector<std::vector<Scalar>>>(dimension.x(), std::vector<std::vector<Scalar>>(dimension.y(), std::vector<Scalar>(dimension.z())));
		m_z[level] = std::vector<std::vector<std::vector<Scalar>>>(dimension.x(), std::vector<std::vector<Scalar>>(dimension.y(), std::vector<Scalar>(dimension.z())));
		m_Adiag[level] = std::vector<std::vector<std::vector<Scalar>>>(dimension.x(), std::vector<std::vector<Scalar>>(dimension.y(), std::vector<Scalar>(dimension.z())));
		m_Acoef[level] = std::vector<std::vector<std::vector<VectorX<6>>>>(dimension.x(), std::vector<std::vector<VectorX<6>>>(dimension.y(), std::vector<VectorX<6>>(dimension.z())));
	}

	m_oldRZ = 0.0;
	m_newRZ = 0.0;
	m_alpha = 0.0;
	m_beta = 0.0;
}

Scalar CGSolver<3>::Reduce(const std::vector<std::vector<std::vector<Scalar>>>& vec0, const std::vector<std::vector<std::vector<Scalar>>>& vec1) const {
	Scalar sum = 0.0;
#pragma omp parallel for reduction(+:sum)
	for (Integer xi = 0; xi < vec0.size(); ++xi) {
		for (Integer yi = 0; yi < vec0[xi].size(); ++yi) {
			for (Integer zi = 0; zi < vec0[xi][yi].size(); ++zi) {
				sum += vec0[xi][yi][zi] * vec1[xi][yi][zi];
			}
		}
	}
	return sum;
}

Scalar CGSolver<3>::GetApproxLiquidPhi(Integer xi, Integer yi, Integer zi, Integer level) const {
	if (level < 0) throw std::invalid_argument("Negative level is not accepted");
	if (level == 0) return m_grid.liquidPhi[xi][yi][zi];
	Integer factor = 1 << level;
	Integer baseX = factor * xi + (1 << (level - 1)) - 1;
	Integer baseY = factor * yi + (1 << (level - 1)) - 1;
	Integer baseZ = factor * zi + (1 << (level - 1)) - 1;
	return MathUtils::Trilerp(m_grid.liquidPhi[baseX][baseY][baseZ], m_grid.liquidPhi[baseX + 1][baseY][baseZ],
						   	  m_grid.liquidPhi[baseX][baseY + 1][baseZ], m_grid.liquidPhi[baseX + 1][baseY + 1][baseZ],
							  m_grid.liquidPhi[baseX][baseY][baseZ + 1], m_grid.liquidPhi[baseX + 1][baseY][baseZ + 1],
							  m_grid.liquidPhi[baseX][baseY + 1][baseZ + 1], m_grid.liquidPhi[baseX + 1][baseY + 1][baseZ + 1], 0.5, 0.5, 0.5);
}

Scalar CGSolver<3>::GetApproxWx(Integer xi, Integer yi, Integer zi, Integer level) const {
	if (level < 0) throw std::invalid_argument("Negative level is not accepted");
	if (level == 0) return m_grid.liquidWx[xi][yi][zi];
	Integer factor = 1 << level;
	Integer baseX = factor * xi;
	Integer baseY = factor * yi + (1 << (level - 1)) - 1;
	Integer baseZ = factor * zi + (1 << (level - 1)) - 1;
	return MathUtils::Bilerp(m_grid.liquidWx[baseX][baseY][baseZ], m_grid.liquidWx[baseX][baseY + 1][baseZ],
							 m_grid.liquidWx[baseX][baseY][baseZ + 1], m_grid.liquidWx[baseX][baseY + 1][baseZ + 1], 0.5, 0.5);
}

Scalar CGSolver<3>::GetApproxWy(Integer xi, Integer yi, Integer zi, Integer level) const {
	if (level < 0) throw std::invalid_argument("Negative level is not accepted");
	if (level == 0) return m_grid.liquidWy[xi][yi][zi];
	Integer factor = 1 << level;
	Integer baseX = factor * xi + (1 << (level - 1)) - 1;
	Integer baseY = factor * yi;
	Integer baseZ = factor * zi + (1 << (level - 1)) - 1;
	return MathUtils::Bilerp(m_grid.liquidWy[baseX][baseY][baseZ], m_grid.liquidWy[baseX + 1][baseY][baseZ],
							 m_grid.liquidWy[baseX][baseY][baseZ + 1], m_grid.liquidWy[baseX + 1][baseY][baseZ + 1], 0.5, 0.5);
}

Scalar CGSolver<3>::GetApproxWz(Integer xi, Integer yi, Integer zi, Integer level) const {
	if (level < 0) throw std::invalid_argument("Negative level is not accepted");
	if (level == 0) return m_grid.liquidWz[xi][yi][zi];
	Integer factor = 1 << level;
	Integer baseX = factor * xi + (1 << (level - 1)) - 1;
	Integer baseY = factor * yi + (1 << (level - 1)) - 1;
	Integer baseZ = factor * zi;
	return MathUtils::Bilerp(m_grid.liquidWz[baseX][baseY][baseZ], m_grid.liquidWz[baseX + 1][baseY][baseZ],
							 m_grid.liquidWz[baseX][baseY + 1][baseZ], m_grid.liquidWz[baseX + 1][baseY + 1][baseZ], 0.5, 0.5);
}

void CGSolver<3>::BuildLHS(const Scalar& liquidDensity, const Scalar& dt) {
	// build levels
	for (Integer level = 0; level < m_maxLevel; ++level) {
		Integer factor = 1 << level;
		Vector3i dimension = m_grid.dimension / factor;
		Vector3 dx = factor * m_grid.dx;

#pragma omp parallel for
		for (Integer xi = 0; xi < dimension.x(); ++xi) {
			for (Integer yi = 0; yi < dimension.y(); ++yi) {
				for (Integer zi = 0; zi < dimension.z(); ++zi) {
					Scalar centerPhi = GetApproxLiquidPhi(xi, yi, zi, level);
					if (centerPhi > 0.0) continue;

					// right
					if (xi + 1 < dimension.x()) {
						Scalar rightWx = GetApproxWx(xi + 1, yi, zi, level);
						if (rightWx > 0.0) {
							Scalar term = rightWx * dt / (liquidDensity * dx.x() * dx.x());
							Scalar phi = GetApproxLiquidPhi(xi + 1, yi, zi, level);
							if (phi < 0.0) {
								m_Acoef[level][xi][yi][zi][0] -= term;
								m_Adiag[level][xi][yi][zi] += term;
							}
							else {
								Scalar frac = max(0.1, MathUtils::FractionInRegion(centerPhi, phi));
								m_Adiag[level][xi][yi][zi] += term / frac;
							}
						}
					}

					// left
					if (xi - 1 >= 0) {
						Scalar leftWx = GetApproxWx(xi, yi, zi, level);
						if (leftWx > 0.0) {
							Scalar term = leftWx * dt / (liquidDensity * dx.x() * dx.x());
							Scalar phi = GetApproxLiquidPhi(xi - 1, yi, zi, level);
							if (phi < 0.0) {
								m_Acoef[level][xi][yi][zi][1] -= term;
								m_Adiag[level][xi][yi][zi] += term;
							}
							else {
								Scalar frac = max(0.1, MathUtils::FractionInRegion(centerPhi, phi));
								m_Adiag[level][xi][yi][zi] += term / frac;
							}
						}
					}

					// top
					if (yi + 1 < dimension.y()) {
						Scalar topWy = GetApproxWy(xi, yi + 1, zi, level);
						if (topWy > 0.0) {
							Scalar term = topWy * dt / (liquidDensity * dx.y() * dx.y());
							Scalar phi = GetApproxLiquidPhi(xi, yi + 1, zi, level);
							if (phi < 0.0) {
								m_Acoef[level][xi][yi][zi][2] -= term;
								m_Adiag[level][xi][yi][zi] += term;
							}
							else {
								Scalar frac = max(0.1, MathUtils::FractionInRegion(centerPhi, phi));
								m_Adiag[level][xi][yi][zi] += term / frac;
							}
						}
					}

					// bottom
					if (yi - 1 >= 0) {
						Scalar botWy = GetApproxWy(xi, yi, zi, level);
						if (botWy > 0.0) {
							Scalar term = botWy * dt / (liquidDensity * dx.y() * dx.y());
							Scalar phi = GetApproxLiquidPhi(xi, yi - 1, zi, level);
							if (phi < 0.0) {
								m_Acoef[level][xi][yi][zi][3] -= term;
								m_Adiag[level][xi][yi][zi] += term;
							}
							else {
								Scalar frac = max(0.1, MathUtils::FractionInRegion(centerPhi, phi));
								m_Adiag[level][xi][yi][zi] += term / frac;
							}
						}
					}

					// front
					if (zi + 1 < dimension.z()) {
						Scalar frontWz = GetApproxWz(xi, yi, zi + 1, level);
						if (frontWz > 0.0) {
							Scalar term = frontWz * dt / (liquidDensity * dx.z() * dx.z());
							Scalar phi = GetApproxLiquidPhi(xi, yi, zi + 1, level);
							if (phi < 0.0) {
								m_Acoef[level][xi][yi][zi][4] -= term;
								m_Adiag[level][xi][yi][zi] += term;
							}
							else {
								Scalar frac = max(0.1, MathUtils::FractionInRegion(centerPhi, phi));
								m_Adiag[level][xi][yi][zi] += term / frac;
							}
						}
					}

					// back
					if (zi - 1 >= 0) {
						Scalar backWz = GetApproxWz(xi, yi, zi, level);
						if (backWz > 0.0) {
							Scalar term = backWz * dt / (liquidDensity * dx.z() * dx.z());
							Scalar phi = GetApproxLiquidPhi(xi, yi, zi - 1, level);
							if (phi < 0.0) {
								m_Acoef[level][xi][yi][zi][5] -= term;
								m_Adiag[level][xi][yi][zi] += term;
							}
							else {
								Scalar frac = max(0.1, MathUtils::FractionInRegion(centerPhi, phi));
								m_Adiag[level][xi][yi][zi] += term / frac;
							}
						}
					}
				}
			}
		}

	}
}

void CGSolver<3>::Init() {
	m_oldRZ = m_newRZ = m_alpha = m_beta = 0.0;

	for (Integer xi = 0; xi < m_grid.dimension.x(); ++xi) {
		for (Integer yi = 0; yi < m_grid.dimension.y(); ++yi) {
			m_p[xi][yi].assign(m_grid.dimension.z(), Scalar(0.0));
			m_x[xi][yi].assign(m_grid.dimension.z(), Scalar(0.0));
			m_rhs[xi][yi].assign(m_grid.dimension.z(), Scalar(0.0));
		}
	}

#pragma omp parallel for
	for (Integer level = 0; level < m_maxLevel; ++level) {
		Vector3i dimension = m_grid.dimension / (1 << level);

		for (Integer xi = 0; xi < dimension.x(); ++xi) {
			for (Integer yi = 0; yi < dimension.y(); ++yi) {
				m_r[level][xi][yi].assign(dimension.z(), 0.0);
				m_z[level][xi][yi].assign(dimension.z(), 0.0);
				m_Adiag[level][xi][yi].assign(dimension.z(), 0.0);
				m_Acoef[level][xi][yi].assign(dimension.z(), VectorX<6>::Zero());
			}
		}
	}
}

void CGSolver<3>::BuildRHS() {
#pragma omp parallel for
	for (Integer xi = 0; xi < m_grid.dimension.x(); ++xi) {
		for (Integer yi = 0; yi < m_grid.dimension.y(); ++yi) {
			for (Integer zi = 0; zi < m_grid.dimension.z(); ++zi) {
				m_rhs[xi][yi][zi] = 0.0;
				if (m_grid.liquidPhi[xi][yi][zi] > 0.0) continue;
				m_rhs[xi][yi][zi] -= m_grid.liquidWx[xi + 1][yi][zi] * m_grid.vx[xi + 1][yi][zi] / m_grid.dx.x();
				m_rhs[xi][yi][zi] += m_grid.liquidWx[xi][yi][zi] * m_grid.vx[xi][yi][zi] / m_grid.dx.x();
				m_rhs[xi][yi][zi] -= m_grid.liquidWy[xi][yi + 1][zi] * m_grid.vy[xi][yi + 1][zi] / m_grid.dx.y();
				m_rhs[xi][yi][zi] += m_grid.liquidWy[xi][yi][zi] * m_grid.vy[xi][yi][zi] / m_grid.dx.y();
				m_rhs[xi][yi][zi] -= m_grid.liquidWz[xi][yi][zi + 1] * m_grid.vz[xi][yi][zi + 1] / m_grid.dx.z();
				m_rhs[xi][yi][zi] += m_grid.liquidWz[xi][yi][zi] * m_grid.vz[xi][yi][zi] / m_grid.dx.z();
			}
		}
	}
}

Scalar CGSolver<3>::GetAx(const std::vector<std::vector<std::vector<Scalar>>>& grid, Integer xi, Integer yi, Integer zi, Integer level) const {
	Vector3i dimension = m_grid.dimension / (1 << level);
#ifdef _DEBUG
	if (xi < 0 || yi < 0 || zi < 0 || xi >= dimension.x() || yi >= dimension.y() || zi >= dimension.z() || level < 0 || level >= m_maxLevel) {
		throw std::invalid_argument("Get invalid or unexpected argument here");
	}

#endif
	const VectorX<6>& coef = m_Acoef[level][xi][yi][zi];
	Scalar result = grid[xi][yi][zi] * m_Adiag[level][xi][yi][zi];
	if (xi + 1 < dimension.x()) result += grid[xi + 1][yi][zi] * coef[0];
	if (xi - 1 >= 0) result += grid[xi - 1][yi][zi] * coef[1];
	if (yi + 1 < dimension.y()) result += grid[xi][yi + 1][zi] * coef[2];
	if (yi - 1 >= 0) result += grid[xi][yi - 1][zi] * coef[3];
	if (zi + 1 < dimension.z()) result += grid[xi][yi][zi + 1] * coef[4];
	if (zi - 1 >= 0) result += grid[xi][yi][zi - 1] * coef[5];
	return result;
}

void CGSolver<3>::ApplyPreconditioner() {
	switch (m_type) {
	case PrecondType::None: {
		for (Integer xi = 0; xi < m_grid.dimension.x(); ++xi) {
			for (Integer yi = 0; yi < m_grid.dimension.y(); ++yi) {
				for (Integer zi = 0; zi < m_grid.dimension.z(); ++zi) {
					if (m_grid.liquidPhi[xi][yi][zi] > 0.0) continue;
					m_z[0][xi][yi][zi] = m_r[0][xi][yi][zi];
				}
			}
		}
		break;
	}
	case PrecondType::Jacobi: {
		for (Integer xi = 0; xi < m_grid.dimension.x(); ++xi) {
			for (Integer yi = 0; yi < m_grid.dimension.y(); ++yi) {
				for (Integer zi = 0; zi < m_grid.dimension.z(); ++zi) {
					const Scalar& denom = m_Adiag[0][xi][yi][zi];
					if (MathUtils::IsSmall(denom)) continue;
					m_z[0][xi][yi][zi] = m_r[0][xi][yi][zi] / denom;
				}
			}
		}
		break;
	}
	case PrecondType::Multigrid: {
		for (Integer level = 1; level < m_maxLevel; ++level) {
			for (Integer iter = 0; iter < m_numPreSmooth; ++iter) {
				Smooth(level - 1, 0);
				Smooth(level - 1, 1);
			}

			Restrict(level);
		}

		for (Integer iter = 0; iter < m_numFinalSmooth; ++iter) {
			Smooth(m_maxLevel - 1, 0);
			Smooth(m_maxLevel - 1, 1);
		}

		for (Integer level = m_maxLevel - 2; level >= 0; --level) {
			Prolong(level);
			for (Integer iter = 0; iter < m_numPostSmooth; ++iter) {
				Smooth(level, 0);
				Smooth(level, 1);
			}
		}
		break;
	}
	default: {
		throw std::runtime_error("Failed to recognize the preconditioner type");
	}
	}
}

void  CGSolver<3>::Prolong(Integer toLevel) {
	Integer fromLevel = toLevel + 1;
	Vector3i dimension = m_grid.dimension / (1 << fromLevel);

#pragma omp parallel for 
	for (Integer xi = 0; xi < dimension.x(); ++xi) {
		for (Integer yi = 0; yi < dimension.y(); ++yi) {
			for (Integer zi = 0; zi < dimension.z(); ++zi) {
				Integer baseX = xi << 1, baseY = yi << 1, baseZ = zi << 1;
				m_z[toLevel][baseX][baseY][baseZ] += m_z[fromLevel][xi][yi][zi];
				m_z[toLevel][baseX][baseY][baseZ + 1] += m_z[fromLevel][xi][yi][zi];
				m_z[toLevel][baseX][baseY + 1][baseZ] += m_z[fromLevel][xi][yi][zi];
				m_z[toLevel][baseX][baseY + 1][baseZ + 1] += m_z[fromLevel][xi][yi][zi];
				m_z[toLevel][baseX + 1][baseY][baseZ] += m_z[fromLevel][xi][yi][zi];
				m_z[toLevel][baseX + 1][baseY][baseZ + 1] += m_z[fromLevel][xi][yi][zi];
				m_z[toLevel][baseX + 1][baseY + 1][baseZ] += m_z[fromLevel][xi][yi][zi];
				m_z[toLevel][baseX + 1][baseY + 1][baseZ + 1] += m_z[fromLevel][xi][yi][zi];
			}
		}
	}
}

void CGSolver<3>::Restrict(Integer toLevel) {
	Integer fromLevel = toLevel - 1;
	Vector3i dimension = m_grid.dimension / (1 << toLevel);
#pragma omp parallel for 
	for (Integer xi = 0; xi < dimension.x(); ++xi) {
		for (Integer yi = 0; yi < dimension.y(); ++yi) {
			for (Integer zi = 0; zi < dimension.z(); ++zi) {
				m_r[toLevel][xi][yi][zi] = 0.0;
				m_z[toLevel][xi][yi][zi] = 0.0;

				Integer baseX = xi << 1, baseY = yi << 1, baseZ = zi << 1;
				m_r[toLevel][xi][yi][zi] += 0.125 * (m_r[fromLevel][baseX][baseY][baseZ] - GetAx(m_z[fromLevel], baseX, baseY, baseZ, fromLevel));
				m_r[toLevel][xi][yi][zi] += 0.125 * (m_r[fromLevel][baseX][baseY][baseZ + 1] - GetAx(m_z[fromLevel], baseX, baseY, baseZ + 1, fromLevel));
				m_r[toLevel][xi][yi][zi] += 0.125 * (m_r[fromLevel][baseX][baseY + 1][baseZ] - GetAx(m_z[fromLevel], baseX, baseY + 1, baseZ, fromLevel));
				m_r[toLevel][xi][yi][zi] += 0.125 * (m_r[fromLevel][baseX][baseY + 1][baseZ + 1] - GetAx(m_z[fromLevel], baseX, baseY + 1, baseZ + 1, fromLevel));
				m_r[toLevel][xi][yi][zi] += 0.125 * (m_r[fromLevel][baseX + 1][baseY][baseZ] - GetAx(m_z[fromLevel], baseX + 1, baseY, baseZ, fromLevel));
				m_r[toLevel][xi][yi][zi] += 0.125 * (m_r[fromLevel][baseX + 1][baseY][baseZ + 1] - GetAx(m_z[fromLevel], baseX + 1, baseY, baseZ + 1, fromLevel));
				m_r[toLevel][xi][yi][zi] += 0.125 * (m_r[fromLevel][baseX + 1][baseY + 1][baseZ] - GetAx(m_z[fromLevel], baseX + 1, baseY + 1, baseZ, fromLevel));
				m_r[toLevel][xi][yi][zi] += 0.125 * (m_r[fromLevel][baseX + 1][baseY + 1][baseZ + 1] - GetAx(m_z[fromLevel], baseX + 1, baseY + 1, baseZ + 1, fromLevel));
			}
		}
	}
}

void CGSolver<3>::Smooth(Integer level, Integer phase) {
	Vector3i dimension = m_grid.dimension / (1 << level);

#pragma omp parallel for
	for (Integer xi = 0; xi < dimension.x(); ++xi) {
		for (Integer yi = 0; yi < dimension.y(); ++yi) {
			for (Integer zi = 0; zi < dimension.z(); ++zi) {
				if (((xi + yi + zi) & 0b1) == phase) {
					if (!MathUtils::IsSmall(m_Adiag[level][xi][yi][zi])) {
						Scalar term = (m_r[level][xi][yi][zi] - GetAx(m_z[level], xi, yi, zi, level)) / m_Adiag[level][xi][yi][zi];
						m_z[level][xi][yi][zi] = 0.35 * m_z[level][xi][yi][zi] + 0.65 * term;
					}
				}
			}
		}
	}
}

void CGSolver<3>::Solve(const Scalar& liquidDensity, const Scalar& dt) {
	// Initialize
	Init();
	BuildRHS();
	BuildLHS(liquidDensity, dt);

	for (Integer xi = 0; xi < m_grid.dimension.x(); ++xi) {
		for (Integer yi = 0; yi < m_grid.dimension.y(); ++yi) {
			for (Integer zi = 0; zi < m_grid.dimension.z(); ++zi) {
				if (m_grid.liquidPhi[xi][yi][zi] > 0.0) continue;
				m_r[0][xi][yi][zi] = m_rhs[xi][yi][zi];
			}
		}
	}
	ApplyPreconditioner();
	for (Integer xi = 0; xi < m_grid.dimension.x(); ++xi) {
		for (Integer yi = 0; yi < m_grid.dimension.y(); ++yi) {
			for (Integer zi = 0; zi < m_grid.dimension.z(); ++zi) {
				if (m_grid.liquidPhi[xi][yi][zi] > 0.0) continue;
				m_p[xi][yi][zi] = m_z[0][xi][yi][zi];
			}
		}
	}

	m_oldRZ = Reduce(m_r[0], m_z[0]);
	Scalar residual = std::sqrt(Reduce(m_r[0], m_r[0]));

	// start real solving
	for (Integer iter = 0; iter < m_maxIteration; ++iter) {
		if (residual < m_tolerance) break;

		Scalar denom = 0.0;
		for (Integer xi = 1; xi < m_grid.dimension.x() - 1; ++xi) {
			for (Integer yi = 1; yi < m_grid.dimension.y() - 1; ++yi) {
				for (Integer zi = 1; zi < m_grid.dimension.z() - 1; ++zi) {
					if (m_grid.liquidPhi[xi][yi][zi] > 0.0) continue;
					denom += m_p[xi][yi][zi] * GetAx(m_p, xi, yi, zi, 0);
				}
			}
		}
		m_alpha = m_oldRZ / denom;
#ifdef _DEBUG
		ASSERT_MSG(!std::isnan(m_alpha) && !std::isinf(m_alpha),
			"Receive alpha=" + std::to_string(m_alpha) + " with oldRZ=" + std::to_string(m_oldRZ) + " and denom=" + std::to_string(denom));
#endif

		residual = 0.0;
#pragma omp parallel for  schedule(dynamic)
		for (Integer xi = 1; xi < m_grid.dimension.x() - 1; ++xi) {
			for (Integer yi = 1; yi < m_grid.dimension.y() - 1; ++yi) {
				for (Integer zi = 1; zi < m_grid.dimension.z() - 1; ++zi) {
					if (m_grid.liquidPhi[xi][yi][zi] > 0.0) continue;

					m_x[xi][yi][zi] += m_alpha * m_p[xi][yi][zi];
					m_r[0][xi][yi][zi] -= m_alpha * GetAx(m_p, xi, yi, zi, 0);
					residual += m_r[0][xi][yi][zi] * m_r[0][xi][yi][zi];
				}
			}
		}
		residual = std::sqrt(residual);

		ApplyPreconditioner();
		m_newRZ = Reduce(m_r[0], m_z[0]);
		m_beta = m_newRZ / m_oldRZ;

#ifdef _DEBUG
		ASSERT_MSG(!std::isnan(m_beta) && !std::isinf(m_beta),
			"Receive beta=" + std::to_string(m_beta) + " with oldRZ=" + std::to_string(m_oldRZ) + " and newRZ=" + std::to_string(m_newRZ));
#endif

#pragma omp parallel for  schedule(dynamic)
		for (Integer xi = 0; xi < m_grid.dimension.x(); ++xi) {
			for (Integer yi = 0; yi < m_grid.dimension.y(); ++yi) {
				for (Integer zi = 0; zi < m_grid.dimension.z(); ++zi) {
					if (m_grid.liquidPhi[xi][yi][zi] > 0.0) continue;
					m_p[xi][yi][zi] = m_z[0][xi][yi][zi] + m_beta * m_p[xi][yi][zi];
				}
			}
		}
		m_oldRZ = m_newRZ;

#ifdef _DEBUG
		std::cout << "Iteration " << iter << ": residual=" << residual << std::endl;
#endif
	}

#ifdef _DEBUG
	ASSERT_MSG(residual < m_tolerance, "Get too large the residual: " + std::to_string(residual));
#endif

	for (Integer xi = 0; xi < m_grid.dimension.x(); ++xi) {
		for (Integer yi = 0; yi < m_grid.dimension.y(); ++yi) {
			for (Integer zi = 0; zi < m_grid.dimension.z(); ++zi) {
				m_grid.pressure[xi][yi][zi] = m_x[xi][yi][zi];
			}
		}
	}
}

#pragma endregion

#pragma endregion
