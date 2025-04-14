#include "Scene.h"

#include "LiquidModel.h"

template class Scene<2>;
template class Scene<3>;

#pragma region Scene 2D
template <>
Scene<2>::Scene(SceneParameter<2> params)
	: m_params(params), m_model(*this) {
	m_boundary.SetChild(new RectBoundary<2>(Vector2::Zero(), Vector2::Ones() * 0.9, true));
	m_model.SetBoundary(&m_boundary);
	
#ifdef _DEBUG
	ASSERT_MSG(std::abs(m_boundary.Compute(Vector2(0.9, 0.0))) < EPSILON, "Get phi=" + std::to_string(m_boundary.Compute(Vector2(0.9, 0.0))));
	ASSERT_MSG(std::abs(m_boundary.Compute(Vector2(0.0, 0.9))) < EPSILON, "Get phi=" + std::to_string(m_boundary.Compute(Vector2(0.0, 0.9))));
	ASSERT_MSG(m_boundary.Compute(Vector2(0.0, 0.0)) > 0.0, "Get phi=" + std::to_string(m_boundary.Compute(Vector2(0.0, 0.0))));
	ASSERT_MSG(m_boundary.Compute(Vector2(1.0, 0.0)) < 1.0, "Get phi=" + std::to_string(m_boundary.Compute(Vector2(1.0, 0.0))));
#endif

}


#pragma endregion

#pragma region Scene 3D
template <>
Scene<3>::Scene(SceneParameter<3> params)
	: m_params(params), m_model(*this) {
	m_boundary.SetChild(new RectBoundary<3>(Vector3::Zero(), Vector3::Ones() * 0.9, true));
	m_model.SetBoundary(&m_boundary);

#ifdef _DEBUG
	ASSERT_MSG(std::abs(m_boundary.Compute(Vector3(0.9, 0.0, 0.9))) < EPSILON, "Get phi=" + std::to_string(m_boundary.Compute(Vector3(0.9, 0.0, 0.9))));
	ASSERT_MSG(std::abs(m_boundary.Compute(Vector3(0.0, 0.9, 0.0))) < EPSILON, "Get phi=" + std::to_string(m_boundary.Compute(Vector3(0.0, 0.9, 0.0))));
	ASSERT_MSG(m_boundary.Compute(Vector3(0.0, 0.0, 0.0)) > 0.0, "Get phi=" + std::to_string(m_boundary.Compute(Vector3(0.0, 0.0, 0.0))));
	ASSERT_MSG(m_boundary.Compute(Vector3(1.0, 0.0, 0.0)) < 1.0, "Get phi=" + std::to_string(m_boundary.Compute(Vector3(1.0, 0.0, 0.0))));
#endif

}


#pragma endregion