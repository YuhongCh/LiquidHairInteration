#include "HairModel.h"

template class HairModel<3>;

template <>
HairModel<3>::HairModel(const Scene<3>& scene)
	: m_scene(scene), m_params(scene.GetHairParameter()), m_ps(m_params.radius) {
	m_ps.Resize(m_params.numParticles);

	m_hairOffset.resize(m_params.numHairs);
	m_hairState.reserve(m_params.numHairs);
	for (Integer i = 0; i < m_params.numHairs; ++i) {
		m_hairState.emplace_back(*this, i);
	}

	Integer particleIndex = 0;
	for (Integer hairIndex = 0; hairIndex < m_params.numHairs; ++hairIndex) {
		m_hairOffset[hairIndex] = particleIndex;
		for (Integer vertexIndex = 0; vertexIndex < m_params.hairInfo[hairIndex].size(); ++vertexIndex) {
			Particle<3>& particle = m_ps.GetParticle(particleIndex);
			particle.position = m_params.hairInfo[hairIndex][vertexIndex].position;
			particle.velocity = m_params.hairInfo[hairIndex][vertexIndex].velocity;
			particle.isFixed = m_params.hairInfo[hairIndex][vertexIndex].isFixed;
			++particleIndex;
		}
	}

    // initialize render utilities
    m_renderProgram = Renderer::CreateShaderProgram("Shader/Scene3D_vertex.glsl", "Shader/Scene3D_fragment.glsl");
    m_minCoordLoc = glGetUniformLocation(m_renderProgram, "minCoord");
    m_maxCoordLoc = glGetUniformLocation(m_renderProgram, "maxCoord");
    m_pointSizeLoc = glGetUniformLocation(m_renderProgram, "pointSize");
    m_pixelColorLoc = glGetUniformLocation(m_renderProgram, "pixelColor");

    glGenVertexArrays(1, &m_VAO);
    glGenBuffers(1, &m_VBO);
    glGenBuffers(1, &m_EBO);

    std::vector<Integer> indices;
    for (Integer hairIndex = 0; hairIndex < m_params.numHairs; ++hairIndex) {
        Integer start = m_hairOffset[hairIndex];
        Integer end = hairIndex == m_params.numHairs - 1 ? m_params.numParticles : m_hairOffset[hairIndex + 1];
        for (Integer vertexIndex = start; vertexIndex < end - 1; ++vertexIndex) {
            indices.push_back(vertexIndex);
            indices.push_back(vertexIndex + 1);
        }
    }
    m_numRenderLines = indices.size();
    glBindVertexArray(m_VAO);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(Integer), indices.data(), GL_STATIC_DRAW);

    glBindBuffer(GL_ARRAY_BUFFER, m_VBO);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Particle<3>), (void*)offsetof(Particle<3>, position));
    glBindVertexArray(0);
}

template <>
HairModel<3>::~HairModel() {
    glDeleteBuffers(1, &m_VBO);
    glDeleteVertexArrays(1, &m_VAO);
    glDeleteBuffers(1, &m_EBO);
    glDeleteProgram(m_renderProgram);
}

template <>
void HairModel<3>::Clear() {
	for (Integer hairIndex = 0; hairIndex < m_params.numHairs; ++hairIndex) {
		m_hairState[hairIndex].Clear();
	}
}

template <>
void HairModel<3>::Init() {
	for (Integer hairIndex = 0; hairIndex < m_params.numHairs; ++hairIndex) {
		m_hairState[hairIndex].Compute(true);
	}
}

template <>
void HairModel<3>::PreStep(const Scalar& dt) {
	for (Integer hairIndex = 0; hairIndex < m_params.numHairs; ++hairIndex) {
		m_hairState[hairIndex].Compute(false);
	}
}

template <>
void HairModel<3>::Step(const Scalar& dt) {
	
	
}

template <>
void HairModel<3>::Render() const {
    glUseProgram(m_renderProgram);
    glLineWidth(5.0);

    // Bind and upload the data
    const std::vector<Particle<3>>& particles = m_ps.GetParticles();
    glBindVertexArray(m_VAO);
    glBindBuffer(GL_ARRAY_BUFFER, m_VBO);
    glBufferData(GL_ARRAY_BUFFER, particles.size() * sizeof(Particle<3>), particles.data(), GL_DYNAMIC_DRAW);

    auto [r, g, b, a] = Color::Black().toScalar();
    glUniform3fv(m_minCoordLoc, 1, GetScene().GetMinCoord().data());
    glUniform3fv(m_maxCoordLoc, 1, GetScene().GetMaxCoord().data());
    glUniform1f(m_pointSizeLoc, 5.0);
    glUniform4f(m_pixelColorLoc, r, g, b, a);

    glDrawElements(GL_LINES, m_numRenderLines, GL_UNSIGNED_INT, BUFFER_OFFSET(0));

    // Unbind and clean up the buffers
    glBindVertexArray(0);
}
