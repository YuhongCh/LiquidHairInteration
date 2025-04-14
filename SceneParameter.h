#pragma once

#include "Definition.h"

template <int DIM>
struct SceneParameter {
	Scalar dt;
	Scalar liquidDensity;
	Scalar gravity;

	VectorX<DIM> minCoord;
	VectorX<DIM> maxCoord;
	VectorXi<DIM> gridDimension;
	Integer numParticles;

	static SceneParameter Default() {
		SceneParameter params;
		params.dt = 0.01;
		params.liquidDensity = 1000.0;
		params.gravity = -9.81;
		params.minCoord = -VectorX<DIM>::Ones();
		params.maxCoord = VectorX<DIM>::Ones();
		params.gridDimension = 64 * VectorXi<DIM>::Ones();
		params.numParticles = 5000;
		return params;
	}
};

