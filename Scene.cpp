#include "Scene.h"

#include "LiquidModel.h"
#include "HairModel.h"
#include "Renderer.h"

//template class Scene<2>;
template class Scene<3>;

#pragma region Scene

template <int DIM>
Scene<DIM>::Scene(const std::string& configFile)
	: m_solverParams(configFile), m_sceneParams(configFile), m_liquidParams(configFile), m_hairParams(configFile), 
	m_liquidModel(std::make_unique<LiquidModel<DIM>>(*this)), 
	m_hairModel(std::make_unique<HairModel<DIM>>(*this)) {

	m_boundary.SetChild(new RectBoundary<DIM>(VectorX<DIM>::Zero(), VectorX<DIM>::Ones() * 0.9, true));
	m_liquidModel->SetBoundary(&m_boundary);
}

template <int DIM>
Scene<DIM>::~Scene() = default;

template <int DIM>
void Scene<DIM>::PlaySimulation(const std::string& filename) {
    Renderer& rend = Renderer::GetInstance();
    std::ifstream reader(filename);
    if (!reader.is_open()) {
        throw std::runtime_error("Failed to open file " + filename);
    }

    std::string line;
    std::getline(reader, line);
    m_liquidModel->GetParticleSystem().Resize(std::stoul(line));

    while (std::getline(reader, line)) {
        rend.ClearScreen();
        m_liquidModel->GetParticleSystem().LoadFromLine(line);
        m_liquidModel->Render();
        rend.LoadScreen();
    }

    reader.close();
}

#pragma endregion
