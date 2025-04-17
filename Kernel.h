#pragma once

#include "Definition.h"
#include "MathUtils.h"

struct BSplineKernel {
	static Scalar LinearWeight(Scalar dist) { 
		dist = std::abs(dist);
		return max(0.0f, 1.0f - dist); 
	}

	static Scalar QuadraticWeight(Scalar dist) {
		dist = std::abs(dist);
		if (dist < 0.5f) return 0.75f - dist * dist;
		else if (dist < 1.5f) return 0.5f * (1.5f - dist) * (1.5f - dist);
		return 0.0f;
	}

	static Scalar CubicWeight(Scalar dist) {
		dist = std::abs(dist);
		if (dist < 1.0f) return 0.5f * dist * dist * dist - dist * dist + 2.f / 3.f;
		else if (dist < 2.0f) return (2.f - dist) * (2.f - dist) * (2.f - dist) / 6.f; 
		return 0.0f;
	}

	static Scalar GradQuadraticWeight(Scalar dist) {
		dist = std::abs(dist);
		if (dist < 0.5f) return -2.f * dist;
		else if (dist < 1.5f) return dist - 1.5f;
		return 0.f;
	}

	static Scalar GradLinearWeight(Scalar dist) {
		if (dist > 0.0 && dist < 1.0)  return -1.0f;
		if (dist < 0.0 && dist > -1.0) return 1.0f;
		return 0.0f;
	}

};