#include "LiquidModel.h"
#include "Kernel.h"
#include "MathUtils.h"

template class LiquidModel<2>;
template class LiquidModel<3>;


#pragma region Liquid Model 2D

template <>
LiquidModel<2>::LiquidModel(const Scene<2>& scene)
	: m_scene(scene), m_grid(scene.GetMinCoord(), scene.GetMaxCoord(), scene.GetGridDimension()),
	m_ps(scene.GetNumParticle(), m_grid.dx.x()),
	m_boundary(nullptr), m_ns(m_ps, m_grid),
	m_solver(m_grid, 500, 1e-5, CGSolver<2>::PrecondType::Multigrid, 4, 5, 5, 10){

	m_correctStep = 0;
	m_correctCycle = 8;
	for (Integer i = 0; i < m_ps.NumParticles(); ++i) {
		m_ps.SetRadius(i, m_grid.dx.x());
	}
}

template <>
void LiquidModel<2>::Init() {
	ComputeSolidPhi();

	Vector2 offset = Scalar(0.3) * (m_grid.maxCoord - m_grid.minCoord);
	m_ps.SpawnParticles(m_grid.minCoord + offset, m_grid.maxCoord - offset);
}

void LiquidModel<2>::PreStep(const Scalar& dt) {
	m_grid.Clear();
	ComputeLiquidPhi();
}

template <>
void LiquidModel<2>::Step(const Scalar& dt) {
	Particle2Grid();
	ApplyGravity(dt);
	Projection(dt);
	Extrapolate();
	Grid2Particle();
	AdvectParticle(dt);
	CorrectVolume(dt);
}

template <>
void LiquidModel<2>::PostStep(const Scalar& dt) {
#ifdef _DEBUG
	Scalar div = ComputeDivergence();

	std::cout << "Current Step End: divergence=" << div << std::endl;
#endif
}

template <>
void LiquidModel<2>::ComputeLiquidPhi() {
	#pragma omp parallel for
	for (Integer pi = 0; pi < m_ps.NumParticles(); ++pi) {
		const Vector2& pos = m_ps.GetPosition(pi);
		const Scalar& radius = m_ps.GetRadius(pi);
		Vector2i coord = m_grid.Point2Index(pos);

		Integer deltaX = static_cast<Integer>(std::ceil(radius / m_grid.dx.x()));
		Integer deltaY = static_cast<Integer>(std::ceil(radius / m_grid.dx.y()));
		for (Integer xi = max(0, coord.x() - deltaX); xi <= min(m_grid.dimension.x() - 1, coord.x() + deltaX); ++xi) {
			for (Integer yi = max(0, coord.y() - deltaY); yi <= min(m_grid.dimension.y() - 1, coord.y() + deltaY); ++yi) {
				Vector2 cellPos = m_grid.Index2Point(xi, yi);
				Scalar phi = (cellPos - pos).norm() - radius;

				#pragma omp critical
				{
					m_grid.liquidPhi[xi][yi] = min(m_grid.liquidPhi[xi][yi], phi);
				}
			}
		}
	}
}

template <>
void LiquidModel<2>::ComputeSolidPhi() {
	if (m_boundary == nullptr) return;

	for (Integer xi = 0; xi < m_grid.dimension.x() + 1; ++xi) {
		for (Integer yi = 0; yi < m_grid.dimension.y() + 1; ++yi) {
			Vector2 cellPos = Vector2(xi, yi).cwiseProduct(m_grid.dx) + m_grid.minCoord;
			m_grid.solidPhi[xi][yi] = m_boundary->Compute(cellPos);
		}
	}

	for (Integer xi = 0; xi < m_grid.dimension.x() + 1; ++xi) {
		for (Integer yi = 0; yi < m_grid.dimension.y(); ++yi) {
			m_grid.liquidWx[xi][yi] = 1.0 - MathUtils::FractionInRegion(m_grid.solidPhi[xi][yi + 1], m_grid.solidPhi[xi][yi]);
			m_grid.liquidWx[xi][yi] = MathUtils::Clamp(m_grid.liquidWx[xi][yi], 0.0, 1.0);
		}
	}

	for (Integer xi = 0; xi < m_grid.dimension.x(); ++xi) {
		for (Integer yi = 0; yi < m_grid.dimension.y() + 1; ++yi) {
			m_grid.liquidWy[xi][yi] = 1.0 - MathUtils::FractionInRegion(m_grid.solidPhi[xi + 1][yi], m_grid.solidPhi[xi][yi]);
			m_grid.liquidWy[xi][yi] = MathUtils::Clamp(m_grid.liquidWy[xi][yi], 0.0, 1.0);
		}
	}
}

template <>
void LiquidModel<2>::Particle2Grid() {
	#pragma omp parallel for
	for (Integer pi = 0; pi < m_ps.NumParticles(); ++pi) {
		const Particle<2>& particle = m_ps.GetParticle(pi);

		Vector2 gridPos = (particle.position - m_grid.minCoord).cwiseQuotient(m_grid.dx);
		Vector2i cellIndex = gridPos.cast<Integer>();
		Integer deltaX = static_cast<Integer>(std::ceil(particle.radius / m_grid.dx.x()));
		Integer deltaY = static_cast<Integer>(std::ceil(particle.radius / m_grid.dx.y()));

		// handle wx and vx grid
		for (Integer dx = -deltaX; dx <= deltaX; ++dx) {
			for (Integer dy = -deltaY; dy <= deltaY; ++dy) {
				Integer xi = cellIndex.x() + dx, yi = cellIndex.y() + dy;
				if (xi < 0 || xi > m_grid.dimension.x() || yi < 0 || yi > m_grid.dimension.y() - 1) continue;
				Vector2 cellCenter(Scalar(xi), yi + Scalar(0.5));
				Scalar weight = BSplineKernel::QuadraticWeight(gridPos.x() - cellCenter.x()) *
								BSplineKernel::QuadraticWeight(gridPos.y() - cellCenter.y());
				Vector2 dpos = particle.position - (cellCenter.cwiseProduct(m_grid.dx) + m_grid.minCoord);

				#pragma omp atomic
				m_grid.wx[xi][yi] += weight;

				#pragma omp atomic
				m_grid.vx[xi][yi] += weight * (particle.velocity[0] + particle.affineMatrix.col(0).dot(dpos));
				
			}
		}
		
		// handle wy and vy grid
		for (Integer dx = -deltaX; dx <= deltaX; ++dx) {
			for (Integer dy = -deltaY; dy <= deltaY; ++dy) {
				Integer xi = cellIndex.x() + dx, yi = cellIndex.y() + dy;
				if (xi < 0 || xi > m_grid.dimension.x() - 1 || yi < 0 || yi > m_grid.dimension.y()) continue;
				Vector2 cellCenter(Scalar(xi) + Scalar(0.5), yi);
				Scalar weight = BSplineKernel::QuadraticWeight(gridPos.x() - cellCenter.x()) *
								BSplineKernel::QuadraticWeight(gridPos.y() - cellCenter.y());
				Vector2 dpos = particle.position - (cellCenter.cwiseProduct(m_grid.dx) + m_grid.minCoord);

				#pragma omp atomic
				m_grid.wy[xi][yi] += weight;

				#pragma omp atomic
				m_grid.vy[xi][yi] += weight * (particle.velocity[1] + particle.affineMatrix.col(1).dot(dpos));
			}
		}
	}

	for (Integer xi = 0; xi < m_grid.dimension.x() + 1; ++xi) {
		for (Integer yi = 0; yi < m_grid.dimension.y(); ++yi) {
			m_grid.vx[xi][yi] = !MathUtils::IsSmall(m_grid.wx[xi][yi]) ? m_grid.vx[xi][yi] / m_grid.wx[xi][yi] : 0.0;
		}
	}

	for (Integer xi = 0; xi < m_grid.dimension.x(); ++xi) {
		for (Integer yi = 0; yi < m_grid.dimension.y() + 1; ++yi) {
			m_grid.vy[xi][yi] = !MathUtils::IsSmall(m_grid.wy[xi][yi]) ? m_grid.vy[xi][yi] / m_grid.wy[xi][yi] : 0.0;
		}
	}
}

template <>
void LiquidModel<2>::Grid2Particle() {
	#pragma omp parallel for
	for (Integer pi = 0; pi < m_ps.NumParticles(); ++pi) {
		const Vector2& pos = m_ps.GetPosition(pi);
		const Scalar& radius = m_ps.GetRadius(pi);

		Vector2 gridPos = (pos - m_grid.minCoord).cwiseQuotient(m_grid.dx);
		Vector2i cellIndex(Integer(gridPos.x()), Integer(gridPos.y()));
		Integer deltaX = static_cast<Integer>(std::ceil(radius / m_grid.dx.x()));
		Integer deltaY = static_cast<Integer>(std::ceil(radius / m_grid.dx.y()));


		Scalar velX = 0.0, velY = 0.0;
		Matrix2 affineMatrix = Matrix2::Zero();

		//handle wx and vx grid
		for (Integer dx = -deltaX; dx <= deltaX; ++dx) {
			for (Integer dy = -deltaY; dy <= deltaY; ++dy) {
				Integer xi = cellIndex.x() + dx, yi = cellIndex.y() + dy;
				if (xi < 0 || xi > m_grid.dimension.x() || yi < 0 || yi >= m_grid.dimension.y()) continue;
				Vector2 cellCenter(Scalar(xi), yi + Scalar(0.5));

				Scalar weightX = BSplineKernel::QuadraticWeight(gridPos.x() - cellCenter.x());
				Scalar weightY = BSplineKernel::QuadraticWeight(gridPos.y() - cellCenter.y());
				Scalar gradWeightX = BSplineKernel::GradQuadraticWeight(gridPos.x() - cellCenter.x());
				Scalar gradWeightY = BSplineKernel::GradQuadraticWeight(gridPos.y() - cellCenter.y());
				Scalar weight = weightX * weightY;
				Vector2 gradWeight(weightY * gradWeightX, weightX * gradWeightY);

				velX += m_grid.vx[xi][yi] * weight;
				affineMatrix.col(0) += m_grid.vx[xi][yi] * gradWeight;
			}
		}

		// handle wy and vy grid
		for (Integer dx = -deltaX; dx <= deltaX; ++dx) {
			for (Integer dy = -deltaY; dy <= deltaY; ++dy) {
				Integer xi = cellIndex.x() + dx, yi = cellIndex.y() + dy;
				if (xi < 0 || xi >= m_grid.dimension.x() || yi < 0 || yi > m_grid.dimension.y()) continue;
				Vector2 cellCenter(xi + Scalar(0.5), Scalar(yi));

				Scalar weightX = BSplineKernel::QuadraticWeight(gridPos.x() - cellCenter.x());
				Scalar weightY = BSplineKernel::QuadraticWeight(gridPos.y() - cellCenter.y());
				Scalar gradWeightX = BSplineKernel::GradQuadraticWeight(gridPos.x() - cellCenter.x());
				Scalar gradWeightY = BSplineKernel::GradQuadraticWeight(gridPos.y() - cellCenter.y());
				Scalar weight = weightX * weightY;
				Vector2 gradWeight(weightY * gradWeightX, weightX * gradWeightY);

				velY += m_grid.vy[xi][yi] * weight;
				affineMatrix.col(1) += m_grid.vy[xi][yi] * gradWeight;
			}
		}

		m_ps.SetVelocity(pi, Vector2(velX, velY));
		m_ps.SetAffineMatrix(pi, affineMatrix);
	}
}

template <>
void LiquidModel<2>::SolvePressure(const Scalar& dt) {
	// solve for pressure
	Scalar liquidDensity = m_scene.GetLiquidDensity();
	m_solver.Solve(liquidDensity, dt);

	// update for velocity
	#pragma omp parallel for
	for (Integer xi = 1; xi < m_grid.dimension.x() - 1; ++xi) {
		for (Integer yi = 1; yi < m_grid.dimension.y() - 1; ++yi) {
			Scalar centerPhi = m_grid.liquidPhi[xi][yi];

			if (m_grid.liquidWx[xi][yi] > 0.0) {
				Scalar phi = m_grid.liquidPhi[xi - 1][yi];
				if (centerPhi < 0.0 || phi < 0.0) {
					Scalar frac = 1.0;
					if (centerPhi >= 0.0 || phi >= 0.0) {
						frac = max(0.1, MathUtils::FractionInRegion(centerPhi, phi));
					}

					Scalar gradPressure = (m_grid.pressure[xi][yi] - m_grid.pressure[xi - 1][yi]) / (m_grid.dx.x() * frac);
					m_grid.vx[xi][yi] -= gradPressure * dt / liquidDensity;
				}
			}
			else {
				m_grid.vx[xi][yi] = 0.0;
			}

			if (m_grid.liquidWy[xi][yi] > 0.0) {
				Scalar phi = m_grid.liquidPhi[xi][yi - 1];
				if (centerPhi < 0.0 || phi < 0.0) {
					Scalar frac = 1.0;
					if (centerPhi >= 0.0 || phi >= 0.0) {
						frac = max(0.1, MathUtils::FractionInRegion(centerPhi, phi));
					}

					Scalar gradPressure = (m_grid.pressure[xi][yi] - m_grid.pressure[xi][yi - 1]) / (m_grid.dx.y() * frac);
					m_grid.vy[xi][yi] -= gradPressure * dt / liquidDensity;
				}
			}
			else {
				m_grid.vy[xi][yi] = 0.0;
			}
		}
	}
}

template <>
void LiquidModel<2>::Projection(const Scalar& dt) {
	SolvePressure(dt);

}

template <>
void LiquidModel<2>::ApplyGravity(const Scalar& dt) {
	#pragma omp parallel for
	for (Integer xi = 1; xi < m_grid.dimension.x() - 1; ++xi) {
		for (Integer yi = 1; yi < m_grid.dimension.y(); ++yi) {
			m_grid.vy[xi][yi] += m_scene.GetGravity() * dt;
		}
	}
}

template <>
void LiquidModel<2>::AdvectParticle(const Scalar& dt) {
	#pragma omp parallel for
	for (Integer pi = 0; pi < m_ps.NumParticles(); ++pi) {
		Vector2 vel = m_ps.GetVelocity(pi);
		Vector2 pos = m_ps.GetPosition(pi);

		pos += vel * dt;
		Scalar solidPhi = m_grid.GetSolidPhi(pos);
		if (solidPhi < 0.0) {
			Vector2 grad = m_grid.GetSolidGradient(pos);
			pos -= (solidPhi + EPSILON) * grad;

			vel[0] = (!MathUtils::IsSmall(grad.x())) ? -vel[0] : vel[0];
			vel[1] = (!MathUtils::IsSmall(grad.y())) ? -vel[1] : vel[1];
		}

		m_ps.SetPosition(pi, pos);
	}
	m_ns.ReCompute();
}

template <>
void LiquidModel<2>::Extrapolate() {

	// process x axis extrapolate
	for (Integer xi = 1; xi < m_grid.dimension.x(); ++xi) {
		for (Integer yi = 1; yi < m_grid.dimension.y() - 1; ++yi) {
			m_grid.isValid[xi][yi] = (m_grid.liquidWx[xi][yi] > 0.0);
		}
	}

	for (Integer iter = 0; iter < 5; ++iter) {
		#pragma omp parallel for
		for (Integer xi = 1; xi < m_grid.dimension.x(); ++xi) {
			for (Integer yi = 1; yi < m_grid.dimension.y() - 1; ++yi) {
				m_grid.tempVx[xi][yi] = m_grid.vx[xi][yi];
				if (m_grid.isValid[xi][yi]) continue;

				Scalar sum = 0.0;
				Integer count = 0;
				for (Integer dx = -1; dx <= 1; ++dx) {
					for (Integer dy = -1; dy <= 1; ++dy) {
						if (m_grid.isValid[xi + dx][yi + dy]) {
							sum += m_grid.vx[xi + dx][yi + dy];
							++count;
						}
					}
				}

				if (count > 0) {
					m_grid.tempVx[xi][yi] = sum / count;
					m_grid.isValid[xi][yi] = true;
				}
			}
		}
		std::swap(m_grid.tempVx, m_grid.vx);
	}

	// process x axis extrapolate
	for (Integer xi = 1; xi < m_grid.dimension.x() - 1; ++xi) {
		for (Integer yi = 1; yi < m_grid.dimension.y(); ++yi) {
			m_grid.isValid[xi][yi] = (m_grid.liquidWy[xi][yi] > 0.0);
		}
	}

	for (Integer iter = 0; iter < 5; ++iter) {
		#pragma omp parallel for
		for (Integer xi = 1; xi < m_grid.dimension.x() - 1; ++xi) {
			for (Integer yi = 1; yi < m_grid.dimension.y(); ++yi) {
				m_grid.tempVy[xi][yi] = m_grid.vy[xi][yi];
				if (m_grid.isValid[xi][yi]) continue;

				Scalar sum = 0.0;
				Integer count = 0;

				for (Integer dx = -1; dx <= 1; ++dx) {
					for (Integer dy = -1; dy <= 1; ++dy) {
						if (m_grid.isValid[xi + dx][yi + dy]) {
							sum += m_grid.vy[xi + dx][yi + dy];
							++count;
						}
					}
				}

				if (count > 0) {
					m_grid.tempVy[xi][yi] = sum / count;
					m_grid.isValid[xi][yi] = true;
				}
			}
		}
		std::swap(m_grid.tempVy, m_grid.vy);
	}
}

template <>
Scalar LiquidModel<2>::ComputeDivergence() const {
	Scalar divergence = 0.0;
	Integer count = 0;
	for (Integer xi = 0; xi < m_grid.dimension.x(); ++xi) {
		for (Integer yi = 0; yi < m_grid.dimension.y(); ++yi) {
			if (m_grid.liquidPhi[xi][yi] > 0.0) continue;
			++count;
			Scalar div = m_grid.liquidWx[xi + 1][yi] * m_grid.vx[xi + 1][yi] - m_grid.liquidWx[xi][yi] * m_grid.vx[xi][yi];
			div += m_grid.liquidWy[xi][yi + 1] * m_grid.vy[xi][yi + 1] - m_grid.liquidWy[xi][yi] * m_grid.vy[xi][yi];
			div /= m_grid.dx.x();
			divergence += div;
		}
	}
	return divergence / count;
}

template <>
void LiquidModel<2>::CorrectVolume(Scalar dt) {
	#pragma omp parallel for schedule(dynamic)
	for (Integer pi = 0; pi < m_ps.NumParticles(); ++pi) {
		if (pi % m_correctCycle != m_correctStep) continue;

		Vector2 repel = Vector2::Zero();
		const Vector2& pos = m_ps.GetPosition(pi);

		m_ns.ForEachNeighborParticles(pi, [&](Integer pi, Integer npi) {
			Scalar radius = m_ps.GetRadius(pi);
			Vector2 dir = pos - m_ps.GetPosition(npi);
			Scalar dist = dir.norm();
			Scalar weight = std::sqrt(m_grid.dx.x() * m_grid.dx.y()) * BSplineKernel::QuadraticWeight(dist / radius);
			if (dist > 1e-4 * radius) {
				repel += weight * radius * (dir / dist);
			}
			else {
				repel[0] += radius / dt * (std::rand() & 0xFF) / 255.0;
				repel[1] += radius / dt * (std::rand() & 0xFF) / 255.0;
			}
		});
		Vector2 newPos = pos + repel * dt;
		Scalar solidPhi = m_grid.GetSolidPhi(newPos);
		if (solidPhi < 0.0) {
			Vector2 grad = m_grid.GetSolidGradient(newPos);
			newPos -= solidPhi * grad;
		}
		m_ps.SetBuffer(pi, newPos);
	}

	for (Integer pi = 0; pi < m_ps.NumParticles(); ++pi) {
		if (pi % m_correctCycle != m_correctStep) continue;
		m_ps.SetPosition(pi, m_ps.GetBuffer(pi));
	}

	m_correctStep = (m_correctStep + 1) % m_correctCycle;
	m_ns.ReCompute();
}



#pragma endregion

#pragma region Liquid Model 3D

template <>
LiquidModel<3>::LiquidModel(const Scene<3>& scene)
	: m_scene(scene), m_grid(scene.GetMinCoord(), scene.GetMaxCoord(), scene.GetGridDimension()),
	m_ps(scene.GetNumParticle(), m_grid.dx.x()),
	m_boundary(nullptr), m_ns(m_ps, m_grid),
	m_solver(m_grid, 500, 1e-5, CGSolver<3>::PrecondType::Multigrid, 4) {

	m_correctStep = 0;
	m_correctCycle = 8;
}

template <>
void LiquidModel<3>::Init() {
	ComputeSolidPhi();

	Vector3 offset = Scalar(0.3) * (m_grid.maxCoord - m_grid.minCoord);
	m_ps.SpawnParticles(m_grid.minCoord + offset, m_grid.maxCoord - offset);
}

void LiquidModel<3>::PreStep(const Scalar& dt) {
	m_grid.Clear();
	ComputeLiquidPhi();
}

template <>
void LiquidModel<3>::Step(const Scalar& dt) {
	Particle2Grid();
	ApplyGravity(dt);
	Projection(dt);
	Extrapolate();
	Grid2Particle();
	AdvectParticle(dt);
	CorrectVolume(dt);
}

template <>
void LiquidModel<3>::PostStep(const Scalar& dt) {
#ifdef _DEBUG
	Scalar div = ComputeDivergence();

	std::cout << "Current Step End: divergence=" << div << std::endl;
#endif
}

template <>
void LiquidModel<3>::ComputeLiquidPhi() {
#pragma omp parallel for
	for (Integer pi = 0; pi < m_ps.NumParticles(); ++pi) {
		const Vector3& pos = m_ps.GetPosition(pi);
		const Scalar& radius = m_ps.GetRadius(pi);
		Vector3i coord = m_grid.Point2Index(pos);

		Integer deltaX = static_cast<Integer>(std::ceil(radius / m_grid.dx.x()));
		Integer deltaY = static_cast<Integer>(std::ceil(radius / m_grid.dx.y()));
		Integer deltaZ = static_cast<Integer>(std::ceil(radius / m_grid.dx.z()));
		for (Integer xi = max(0, coord.x() - deltaX); xi <= min(m_grid.dimension.x() - 1, coord.x() + deltaX); ++xi) {
			for (Integer yi = max(0, coord.y() - deltaY); yi <= min(m_grid.dimension.y() - 1, coord.y() + deltaY); ++yi) {
				for (Integer zi = max(0, coord.z() - deltaZ); zi <= min(m_grid.dimension.z() - 1, coord.z() + deltaZ); ++zi) {
					Vector3 cellPos = m_grid.Index2Point(xi, yi, zi);
					Scalar phi = (cellPos - pos).norm() - radius;

					#pragma omp critical
					{
						m_grid.liquidPhi[xi][yi][zi] = min(m_grid.liquidPhi[xi][yi][zi], phi);
					}
				}
			}
		}
	}
}

template <>
void LiquidModel<3>::ComputeSolidPhi() {
	if (m_boundary == nullptr) return;

	for (Integer xi = 0; xi < m_grid.dimension.x() + 1; ++xi) {
		for (Integer yi = 0; yi < m_grid.dimension.y() + 1; ++yi) {
			for (Integer zi = 0; zi < m_grid.dimension.z() + 1; ++zi) {
				Vector3 cellPos = Vector3(xi, yi, zi).cwiseProduct(m_grid.dx) + m_grid.minCoord;
				m_grid.solidPhi[xi][yi][zi] = m_boundary->Compute(cellPos);
			}
		}
	}

	#pragma omp parallel for
	for (Integer xi = 0; xi < m_grid.dimension.x() + 1; ++xi) {
		for (Integer yi = 0; yi < m_grid.dimension.y(); ++yi) {
			for (Integer zi = 0; zi < m_grid.dimension.z(); ++zi) {
				m_grid.liquidWx[xi][yi][zi] = 1.0 - MathUtils::FractionInRegion(m_grid.solidPhi[xi][yi][zi], m_grid.solidPhi[xi][yi + 1][zi],
																				m_grid.solidPhi[xi][yi][zi + 1], m_grid.solidPhi[xi][yi + 1][zi + 1]);
				m_grid.liquidWx[xi][yi][zi] = MathUtils::Clamp(m_grid.liquidWx[xi][yi][zi], 0.0, 1.0);
			}
		}
	}

	#pragma omp parallel for
	for (Integer xi = 0; xi < m_grid.dimension.x(); ++xi) {
		for (Integer yi = 0; yi < m_grid.dimension.y() + 1; ++yi) {
			for (Integer zi = 0; zi < m_grid.dimension.z(); ++zi) {
				m_grid.liquidWy[xi][yi][zi] = 1.0 - MathUtils::FractionInRegion(m_grid.solidPhi[xi][yi][zi], m_grid.solidPhi[xi + 1][yi][zi],
																				m_grid.solidPhi[xi][yi][zi + 1], m_grid.solidPhi[xi + 1][yi][zi + 1]);
				m_grid.liquidWy[xi][yi][zi] = MathUtils::Clamp(m_grid.liquidWy[xi][yi][zi], 0.0, 1.0);
			}
		}
	}

	#pragma omp parallel for
	for (Integer xi = 0; xi < m_grid.dimension.x(); ++xi) {
		for (Integer yi = 0; yi < m_grid.dimension.y(); ++yi) {
			for (Integer zi = 0; zi < m_grid.dimension.z() + 1; ++zi) {
				m_grid.liquidWz[xi][yi][zi] = 1.0 - MathUtils::FractionInRegion(m_grid.solidPhi[xi][yi][zi], m_grid.solidPhi[xi + 1][yi][zi],
																				m_grid.solidPhi[xi][yi + 1][zi], m_grid.solidPhi[xi + 1][yi + 1][zi]);
				m_grid.liquidWz[xi][yi][zi] = MathUtils::Clamp(m_grid.liquidWz[xi][yi][zi], 0.0, 1.0);
			}
		}
	}
}

template <>
void LiquidModel<3>::Particle2Grid() {
	#pragma omp parallel for 
	for (Integer pi = 0; pi < m_ps.NumParticles(); ++pi) {
		const Vector3& pos = m_ps.GetPosition(pi);
		const Vector3& vel = m_ps.GetVelocity(pi);
		const Matrix3& affineMat = m_ps.GetAffineMatrix(pi);
		const Scalar& radius = m_ps.GetRadius(pi);

		Vector3 gridPos = (pos - m_grid.minCoord).cwiseQuotient(m_grid.dx);
		Vector3i cellIndex = gridPos.cast<Integer>();
		Integer deltaX = static_cast<Integer>(std::ceil(radius / m_grid.dx.x()));
		Integer deltaY = static_cast<Integer>(std::ceil(radius / m_grid.dx.y()));
		Integer deltaZ = static_cast<Integer>(std::ceil(radius / m_grid.dx.z()));

		// handle wx and vx grid
		for (Integer dx = -deltaX; dx <= deltaX; ++dx) {
			for (Integer dy = -deltaY; dy <= deltaY; ++dy) {
				for (Integer dz = -deltaZ; dz <= deltaZ; ++dz) {
					Integer xi = cellIndex.x() + dx, yi = cellIndex.y() + dy, zi = cellIndex.z() + dz;

					if (xi < 0 || xi > m_grid.dimension.x() || yi < 0 || yi > m_grid.dimension.y() - 1 || zi < 0 || zi > m_grid.dimension.z() - 1) continue;
					Vector3 cellCenter(Scalar(xi), yi + Scalar(0.5), zi + Scalar(0.5));
					Scalar weight = BSplineKernel::QuadraticWeight(gridPos.x() - cellCenter.x()) *
									BSplineKernel::QuadraticWeight(gridPos.y() - cellCenter.y()) * 
									BSplineKernel::QuadraticWeight(gridPos.z() - cellCenter.z());
					Vector3 dpos = pos - (cellCenter.cwiseProduct(m_grid.dx) + m_grid.minCoord);

#pragma omp atomic
					m_grid.wx[xi][yi][zi] += weight;

#pragma omp atomic
					m_grid.vx[xi][yi][zi] += weight * (vel[0] + affineMat.col(0).dot(dpos));
				}
			}
		}

		// handle wy and vy grid
		for (Integer dx = -deltaX; dx <= deltaX; ++dx) {
			for (Integer dy = -deltaY; dy <= deltaY; ++dy) {
				for (Integer dz = -deltaZ; dz <= deltaZ; ++dz) {
					Integer xi = cellIndex.x() + dx, yi = cellIndex.y() + dy, zi = cellIndex.z() + dz;

					if (xi < 0 || xi > m_grid.dimension.x() - 1 || yi < 0 || yi > m_grid.dimension.y() || zi < 0 || zi > m_grid.dimension.z() - 1) continue;
					Vector3 cellCenter(Scalar(xi) + Scalar(0.5), yi, zi + Scalar(0.5));
					Scalar weight = BSplineKernel::QuadraticWeight(gridPos.x() - cellCenter.x()) *
									BSplineKernel::QuadraticWeight(gridPos.y() - cellCenter.y()) *
									BSplineKernel::QuadraticWeight(gridPos.z() - cellCenter.z());
					Vector3 dpos = pos - (cellCenter.cwiseProduct(m_grid.dx) + m_grid.minCoord);

#pragma omp atomic
					m_grid.wy[xi][yi][zi] += weight;

#pragma omp atomic
					m_grid.vy[xi][yi][zi] += weight * (vel[1] + affineMat.col(1).dot(dpos));
				}
			}
		}

		// handle wz and vz grid
		for (Integer dx = -deltaX; dx <= deltaX; ++dx) {
			for (Integer dy = -deltaY; dy <= deltaY; ++dy) {
				for (Integer dz = -deltaZ; dz <= deltaZ; ++dz) {
					Integer xi = cellIndex.x() + dx, yi = cellIndex.y() + dy, zi = cellIndex.z() + dz;

					if (xi < 0 || xi > m_grid.dimension.x() - 1 || yi < 0 || yi > m_grid.dimension.y() - 1 || zi < 0 || zi > m_grid.dimension.z()) continue;
					Vector3 cellCenter(Scalar(xi) + Scalar(0.5), yi + Scalar(0.5), zi);
					Scalar weight = BSplineKernel::QuadraticWeight(gridPos.x() - cellCenter.x()) *
									BSplineKernel::QuadraticWeight(gridPos.y() - cellCenter.y()) *
									BSplineKernel::QuadraticWeight(gridPos.z() - cellCenter.z());
					Vector3 dpos = pos - (cellCenter.cwiseProduct(m_grid.dx) + m_grid.minCoord);

#pragma omp atomic
					m_grid.wz[xi][yi][zi] += weight;

#pragma omp atomic
					m_grid.vz[xi][yi][zi] += weight * (vel[2] + affineMat.col(2).dot(dpos));
				}
			}
		}

	}

	for (Integer xi = 0; xi < m_grid.dimension.x() + 1; ++xi) {
		for (Integer yi = 0; yi < m_grid.dimension.y(); ++yi) {
			for (Integer zi = 0; zi < m_grid.dimension.z(); ++zi) {
				m_grid.vx[xi][yi][zi] = !MathUtils::IsSmall(m_grid.wx[xi][yi][zi]) ? m_grid.vx[xi][yi][zi] / m_grid.wx[xi][yi][zi] : 0.0;
			}
		}
	}

	for (Integer xi = 0; xi < m_grid.dimension.x(); ++xi) {
		for (Integer yi = 0; yi < m_grid.dimension.y() + 1; ++yi) {
			for (Integer zi = 0; zi < m_grid.dimension.z(); ++zi) {
				m_grid.vy[xi][yi][zi] = !MathUtils::IsSmall(m_grid.wy[xi][yi][zi]) ? m_grid.vy[xi][yi][zi] / m_grid.wy[xi][yi][zi] : 0.0;
			}
		}
	}

	for (Integer xi = 0; xi < m_grid.dimension.x(); ++xi) {
		for (Integer yi = 0; yi < m_grid.dimension.y(); ++yi) {
			for (Integer zi = 0; zi < m_grid.dimension.z() + 1; ++zi) {
				m_grid.vz[xi][yi][zi] = !MathUtils::IsSmall(m_grid.wz[xi][yi][zi]) ? m_grid.vz[xi][yi][zi] / m_grid.wz[xi][yi][zi] : 0.0;
			}
		}
	}
}

template <>
void LiquidModel<3>::Grid2Particle() {
#pragma omp parallel for
	for (Integer pi = 0; pi < m_ps.NumParticles(); ++pi) {
		const Vector3& pos = m_ps.GetPosition(pi);
		const Scalar& radius = m_ps.GetRadius(pi);

		Vector3 gridPos = (pos - m_grid.minCoord).cwiseQuotient(m_grid.dx);
		Vector3i cellIndex(Integer(gridPos.x()), Integer(gridPos.y()), Integer(gridPos.z()));
		Integer deltaX = static_cast<Integer>(std::ceil(radius / m_grid.dx.x()));
		Integer deltaY = static_cast<Integer>(std::ceil(radius / m_grid.dx.y()));
		Integer deltaZ = static_cast<Integer>(std::ceil(radius / m_grid.dx.z()));


		Scalar velX = 0.0, velY = 0.0, velZ = 0.0;
		Matrix3 affineMatrix = Matrix3::Zero();

		//handle wx and vx grid
		for (Integer dx = -deltaX; dx <= deltaX; ++dx) {
			for (Integer dy = -deltaY; dy <= deltaY; ++dy) {
				for (Integer dz = -deltaZ; dz <= deltaZ; ++dz) {
					Integer xi = cellIndex.x() + dx, yi = cellIndex.y() + dy, zi = cellIndex.z() + dz;
					if (xi < 0 || xi > m_grid.dimension.x() || 
						yi < 0 || yi >= m_grid.dimension.y() ||
						zi < 0 || zi >= m_grid.dimension.z()) continue;
					Vector3 cellCenter(Scalar(xi), yi + Scalar(0.5), zi + Scalar(0.5));

					Scalar weightX = BSplineKernel::QuadraticWeight(gridPos.x() - cellCenter.x());
					Scalar weightY = BSplineKernel::QuadraticWeight(gridPos.y() - cellCenter.y());
					Scalar weightZ = BSplineKernel::QuadraticWeight(gridPos.z() - cellCenter.z());
					Scalar gradWeightX = BSplineKernel::GradQuadraticWeight(gridPos.x() - cellCenter.x());
					Scalar gradWeightY = BSplineKernel::GradQuadraticWeight(gridPos.y() - cellCenter.y());
					Scalar gradWeightZ = BSplineKernel::GradQuadraticWeight(gridPos.z() - cellCenter.z());
					
					Scalar weight = weightX * weightY * weightZ;
					Vector3 gradWeight(gradWeightX * weightY * weightZ, weightX * gradWeightY * weightZ, weightX * weightY * gradWeightZ);

					velX += m_grid.vx[xi][yi][zi] * weight;
					affineMatrix.col(0) += m_grid.vx[xi][yi][zi] * gradWeight;
				}
			}
		}

		// handle wy and vy grid
		for (Integer dx = -deltaX; dx <= deltaX; ++dx) {
			for (Integer dy = -deltaY; dy <= deltaY; ++dy) {
				for (Integer dz = -deltaZ; dz <= deltaZ; ++dz) {
				for (Integer dz = -deltaZ; dz <= deltaZ; ++dz) {
					Integer xi = cellIndex.x() + dx, yi = cellIndex.y() + dy, zi = cellIndex.z() + dz;
					if (xi < 0 || xi >= m_grid.dimension.x() ||
						yi < 0 || yi > m_grid.dimension.y() ||
						zi < 0 || zi >= m_grid.dimension.z()) continue;
					Vector3 cellCenter(Scalar(xi) + Scalar(0.5), yi, zi + Scalar(0.5));

					Scalar weightX = BSplineKernel::QuadraticWeight(gridPos.x() - cellCenter.x());
					Scalar weightY = BSplineKernel::QuadraticWeight(gridPos.y() - cellCenter.y());
					Scalar weightZ = BSplineKernel::QuadraticWeight(gridPos.z() - cellCenter.z());
					Scalar gradWeightX = BSplineKernel::GradQuadraticWeight(gridPos.x() - cellCenter.x());
					Scalar gradWeightY = BSplineKernel::GradQuadraticWeight(gridPos.y() - cellCenter.y());
					Scalar gradWeightZ = BSplineKernel::GradQuadraticWeight(gridPos.z() - cellCenter.z());

					Scalar weight = weightX * weightY * weightZ;
					Vector3 gradWeight(gradWeightX * weightY * weightZ, weightX * gradWeightY * weightZ, weightX * weightY * gradWeightZ);

					velY += m_grid.vy[xi][yi][zi] * weight;
					affineMatrix.col(1) += m_grid.vy[xi][yi][zi] * gradWeight;
				}
			}
		}

		// handle wz and vz grid
		for (Integer dx = -deltaX; dx <= deltaX; ++dx) {
			for (Integer dy = -deltaY; dy <= deltaY; ++dy) {
				for (Integer dz = -deltaZ; dz <= deltaZ; ++dz) {
					Integer xi = cellIndex.x() + dx, yi = cellIndex.y() + dy, zi = cellIndex.z() + dz;
					if (xi < 0 || xi >= m_grid.dimension.x() ||
						yi < 0 || yi >= m_grid.dimension.y() ||
						zi < 0 || zi > m_grid.dimension.z()) continue;
					Vector3 cellCenter(Scalar(xi) + Scalar(0.5), yi + Scalar(0.5), zi);

					Scalar weightX = BSplineKernel::QuadraticWeight(gridPos.x() - cellCenter.x());
					Scalar weightY = BSplineKernel::QuadraticWeight(gridPos.y() - cellCenter.y());
					Scalar weightZ = BSplineKernel::QuadraticWeight(gridPos.z() - cellCenter.z());
					Scalar gradWeightX = BSplineKernel::GradQuadraticWeight(gridPos.x() - cellCenter.x());
					Scalar gradWeightY = BSplineKernel::GradQuadraticWeight(gridPos.y() - cellCenter.y());
					Scalar gradWeightZ = BSplineKernel::GradQuadraticWeight(gridPos.z() - cellCenter.z());

					Scalar weight = weightX * weightY * weightZ;
					Vector3 gradWeight(gradWeightX * weightY * weightZ, weightX * gradWeightY * weightZ, weightX * weightY * gradWeightZ);

					velZ += m_grid.vz[xi][yi][zi] * weight;
					affineMatrix.col(2) += m_grid.vz[xi][yi][zi] * gradWeight;
				}
			}
		}

		m_ps.SetVelocity(pi, Vector3(velX, velY, velZ));
		m_ps.SetAffineMatrix(pi, affineMatrix);
	}
}

template <>
void LiquidModel<3>::SolvePressure(const Scalar& dt) {
	// solve for pressure
	Scalar liquidDensity = m_scene.GetLiquidDensity();
	m_solver.Solve(liquidDensity, dt);

	// update for velocity
#pragma omp parallel for
	for (Integer xi = 1; xi < m_grid.dimension.x() - 1; ++xi) {
		for (Integer yi = 1; yi < m_grid.dimension.y() - 1; ++yi) {
			for (Integer zi = 1; zi < m_grid.dimension.z() - 1; ++zi) {
				Scalar centerPhi = m_grid.liquidPhi[xi][yi][zi];

				if (m_grid.liquidWx[xi][yi][zi] > 0.0) {
					Scalar phi = m_grid.liquidPhi[xi - 1][yi][zi];
					if (centerPhi < 0.0 || phi < 0.0) {
						Scalar frac = 1.0;
						if (centerPhi >= 0.0 || phi >= 0.0) {
							frac = max(0.1, MathUtils::FractionInRegion(centerPhi, phi));
						}

						Scalar gradPressure = (m_grid.pressure[xi][yi][zi] - m_grid.pressure[xi - 1][yi][zi]) / (m_grid.dx.x() * frac);
						m_grid.vx[xi][yi][zi] -= gradPressure * dt / liquidDensity;
					}
				}
				else {
					m_grid.vx[xi][yi][zi] = 0.0;
				}

				if (m_grid.liquidWy[xi][yi][zi] > 0.0) {
					Scalar phi = m_grid.liquidPhi[xi][yi - 1][zi];
					if (centerPhi < 0.0 || phi < 0.0) {
						Scalar frac = 1.0;
						if (centerPhi >= 0.0 || phi >= 0.0) {
							frac = max(0.1, MathUtils::FractionInRegion(centerPhi, phi));
						}

						Scalar gradPressure = (m_grid.pressure[xi][yi][zi] - m_grid.pressure[xi][yi - 1][zi]) / (m_grid.dx.y() * frac);
						m_grid.vy[xi][yi][zi] -= gradPressure * dt / liquidDensity;
					}
				}
				else {
					m_grid.vy[xi][yi][zi] = 0.0;
				}

				if (m_grid.liquidWz[xi][yi][zi] > 0.0) {
					Scalar phi = m_grid.liquidPhi[xi][yi][zi - 1];
					if (centerPhi < 0.0 || phi < 0.0) {
						Scalar frac = 1.0;
						if (centerPhi >= 0.0 || phi >= 0.0) {
							frac = max(0.1, MathUtils::FractionInRegion(centerPhi, phi));
						}

						Scalar gradPressure = (m_grid.pressure[xi][yi][zi] - m_grid.pressure[xi][yi][zi - 1]) / (m_grid.dx.z() * frac);
						m_grid.vz[xi][yi][zi] -= gradPressure * dt / liquidDensity;
					}
				}
				else {
					m_grid.vz[xi][yi][zi] = 0.0;
				}
			}
		}
	}
}

template <>
void LiquidModel<3>::Projection(const Scalar& dt) {
	SolvePressure(dt);

}

template <>
void LiquidModel<3>::ApplyGravity(const Scalar& dt) {
#pragma omp parallel for
	for (Integer xi = 1; xi < m_grid.dimension.x() - 1; ++xi) {
		for (Integer yi = 1; yi < m_grid.dimension.y(); ++yi) {
			for (Integer zi = 1; zi < m_grid.dimension.z() - 1; ++zi) {
				m_grid.vy[xi][yi][zi] += m_scene.GetGravity() * dt;
			}
		}
	}
}

template <>
void LiquidModel<3>::AdvectParticle(const Scalar& dt) {
#pragma omp parallel for
	for (Integer pi = 0; pi < m_ps.NumParticles(); ++pi) {
		Vector3 vel = m_ps.GetVelocity(pi);
		Vector3 pos = m_ps.GetPosition(pi);

		pos += vel * dt;
		Scalar solidPhi = m_grid.GetSolidPhi(pos);
		if (solidPhi < 0.0) {
			Vector3 grad = m_grid.GetSolidGradient(pos);
			pos -= (solidPhi + EPSILON) * grad;

			vel[0] = (!MathUtils::IsSmall(grad.x())) ? -vel[0] : vel[0];
			vel[1] = (!MathUtils::IsSmall(grad.y())) ? -vel[1] : vel[1];
			vel[2] = (!MathUtils::IsSmall(grad.z())) ? -vel[2] : vel[2];
		}

		m_ps.SetPosition(pi, pos);
	}
	m_ns.ReCompute();
}

template <>
void LiquidModel<3>::Extrapolate() {

	// process x axis extrapolate
	for (Integer xi = 1; xi < m_grid.dimension.x(); ++xi) {
		for (Integer yi = 1; yi < m_grid.dimension.y() - 1; ++yi) {
			for (Integer zi = 1; zi < m_grid.dimension.z() - 1; ++zi) {
				m_grid.isValid[xi][yi][zi] = (m_grid.liquidWx[xi][yi][zi] > 0.0);
			}
		}
	}

	for (Integer iter = 0; iter < 5; ++iter) {
#pragma omp parallel for
		for (Integer xi = 1; xi < m_grid.dimension.x(); ++xi) {
			for (Integer yi = 1; yi < m_grid.dimension.y() - 1; ++yi) {
				for (Integer zi = 1; zi < m_grid.dimension.z() - 1; ++zi) {
					m_grid.tempVx[xi][yi][zi] = m_grid.vx[xi][yi][zi];
					if (m_grid.isValid[xi][yi][zi]) continue;

					Scalar sum = 0.0;
					Integer count = 0;
					for (Integer dx = -1; dx <= 1; ++dx) {
						for (Integer dy = -1; dy <= 1; ++dy) {
							for (Integer dz = -1; dz <= 1; ++dz) {
								if (m_grid.isValid[xi + dx][yi + dy][zi + dz]) {
									sum += m_grid.vx[xi + dx][yi + dy][zi + dz];
									++count;
								}
							}
						}
					}

					if (count > 0) {
						m_grid.tempVx[xi][yi][zi] = sum / count;
						m_grid.isValid[xi][yi][zi] = true;
					}
				}
			}
		}
		std::swap(m_grid.tempVx, m_grid.vx);
	}

	// process y axis extrapolate
	for (Integer xi = 1; xi < m_grid.dimension.x() - 1; ++xi) {
		for (Integer yi = 1; yi < m_grid.dimension.y(); ++yi) {
			for (Integer zi = 1; zi < m_grid.dimension.z() - 1; ++zi) {
				m_grid.isValid[xi][yi][zi] = (m_grid.liquidWy[xi][yi][zi] > 0.0);
			}
		}
	}

	for (Integer iter = 0; iter < 5; ++iter) {
#pragma omp parallel for
		for (Integer xi = 1; xi < m_grid.dimension.x() - 1; ++xi) {
			for (Integer yi = 1; yi < m_grid.dimension.y(); ++yi) {
				for (Integer zi = 1; zi < m_grid.dimension.z() - 1; ++zi) {
					m_grid.tempVy[xi][yi][zi] = m_grid.vy[xi][yi][zi];
					if (m_grid.isValid[xi][yi][zi]) continue;

					Scalar sum = 0.0;
					Integer count = 0;
					for (Integer dx = -1; dx <= 1; ++dx) {
						for (Integer dy = -1; dy <= 1; ++dy) {
							for (Integer dz = -1; dz <= 1; ++dz) {
								if (m_grid.isValid[xi + dx][yi + dy][zi + dz]) {
									sum += m_grid.vy[xi + dx][yi + dy][zi + dz];
									++count;
								}
							}
						}
					}

					if (count > 0) {
						m_grid.tempVy[xi][yi][zi] = sum / count;
						m_grid.isValid[xi][yi][zi] = true;
					}
				}
			}
		}
		std::swap(m_grid.tempVy, m_grid.vy);
	}

	// process z axis extrapolate
	for (Integer xi = 1; xi < m_grid.dimension.x() - 1; ++xi) {
		for (Integer yi = 1; yi < m_grid.dimension.y() - 1; ++yi) {
			for (Integer zi = 1; zi < m_grid.dimension.z(); ++zi) {
				m_grid.isValid[xi][yi][zi] = (m_grid.liquidWz[xi][yi][zi] > 0.0);
			}
		}
	}

	for (Integer iter = 0; iter < 5; ++iter) {
#pragma omp parallel for
		for (Integer xi = 1; xi < m_grid.dimension.x() - 1; ++xi) {
			for (Integer yi = 1; yi < m_grid.dimension.y() - 1; ++yi) {
				for (Integer zi = 1; zi < m_grid.dimension.z(); ++zi) {
					m_grid.tempVz[xi][yi][zi] = m_grid.vz[xi][yi][zi];
					if (m_grid.isValid[xi][yi][zi]) continue;

					Scalar sum = 0.0;
					Integer count = 0;
					for (Integer dx = -1; dx <= 1; ++dx) {
						for (Integer dy = -1; dy <= 1; ++dy) {
							for (Integer dz = -1; dz <= 1; ++dz) {
								if (m_grid.isValid[xi + dx][yi + dy][zi + dz]) {
									sum += m_grid.vz[xi + dx][yi + dy][zi + dz];
									++count;
								}
							}
						}
					}

					if (count > 0) {
						m_grid.tempVz[xi][yi][zi] = sum / count;
						m_grid.isValid[xi][yi][zi] = true;
					}
				}
			}
		}
		std::swap(m_grid.tempVz, m_grid.vz);
	}
}

template <>
Scalar LiquidModel<3>::ComputeDivergence() const {
	Scalar divergence = 0.0;
	Integer count = 0;
	for (Integer xi = 0; xi < m_grid.dimension.x(); ++xi) {
		for (Integer yi = 0; yi < m_grid.dimension.y(); ++yi) {
			for (Integer zi = 0; zi < m_grid.dimension.z(); ++zi) {
				if (m_grid.liquidPhi[xi][yi][zi] > 0.0) continue;
				++count;
				Scalar div = m_grid.liquidWx[xi + 1][yi][zi] * m_grid.vx[xi + 1][yi][zi] - m_grid.liquidWx[xi][yi][zi] * m_grid.vx[xi][yi][zi];
				div += m_grid.liquidWy[xi][yi + 1][zi] * m_grid.vy[xi][yi + 1][zi] - m_grid.liquidWy[xi][yi][zi] * m_grid.vy[xi][yi][zi];
				div += m_grid.liquidWz[xi][yi][zi + 1] * m_grid.vz[xi][yi][zi + 1] - m_grid.liquidWz[xi][yi][zi] * m_grid.vz[xi][yi][zi];
				div /= m_grid.dx.x();
				divergence += div;
			}
		}
	}
	return divergence / count;
}

template <>
void LiquidModel<3>::CorrectVolume(Scalar dt) {
	Scalar weightFactor = std::sqrt(m_grid.dx.x() * m_grid.dx.y() * m_grid.dx.z());

#pragma omp parallel for schedule(dynamic)
	for (Integer pi = 0; pi < m_ps.NumParticles(); ++pi) {
		if (pi % m_correctCycle != m_correctStep) continue;

		Vector3 repel = Vector3::Zero();
		const Vector3& pos = m_ps.GetPosition(pi);

		m_ns.ForEachNeighborParticles(pi, [&](Integer pi, Integer npi) {
			Scalar radius = m_ps.GetRadius(pi);
			Vector3 dir = pos - m_ps.GetPosition(npi);
			Scalar dist = dir.norm();
			Scalar weight = weightFactor * BSplineKernel::QuadraticWeight(dist / radius);
			if (dist > 1e-4 * radius) {
				repel += weight * radius * (dir / dist);
			}
			else {
				thread_local std::mt19937 rng(std::random_device{}());
				std::uniform_real_distribution<Scalar> distRand(0.0, 1.0);

				repel[0] += radius / dt * distRand(rng);
				repel[1] += radius / dt * distRand(rng);
				repel[2] += radius / dt * distRand(rng);
			}
			});
		Vector3 newPos = pos + repel * dt;
		Scalar solidPhi = m_grid.GetSolidPhi(newPos);
		if (solidPhi < 0.0) {
			Vector3 grad = m_grid.GetSolidGradient(newPos);
			newPos -= solidPhi * grad;
		}
		m_ps.SetBuffer(pi, newPos);
	}

	for (Integer pi = 0; pi < m_ps.NumParticles(); ++pi) {
		if (pi % m_correctCycle != m_correctStep) continue;
		m_ps.SetPosition(pi, m_ps.GetBuffer(pi));
	}

	m_correctStep = (m_correctStep + 1) % m_correctCycle;
	m_ns.ReCompute();
}



#pragma endregion