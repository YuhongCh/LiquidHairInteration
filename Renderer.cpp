#include "Renderer.h"
#include "ParticleSystem.h"

Renderer::Renderer(Integer screenWidth, Integer screenHeight, const std::string& screenTitle)
    : m_screenWidth(screenWidth), m_screenHeight(screenHeight), m_window(nullptr) {
    if (!glfwInit()) {
        throw std::runtime_error("Failed to initialize GLFW");
    }

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
}

Renderer::~Renderer() {
    if (m_window != nullptr) {
        glfwDestroyWindow(m_window);
    }
    glfwTerminate();
}

Renderer& Renderer::GetInstance() {
    static Renderer renderer;
    return renderer;
}

void Renderer::InstantiateWindow(const std::string& title, const Integer& width, const Integer& height) {
    if (m_window) {
        glfwDestroyWindow(m_window);
        m_window = nullptr;
    }

    m_screenWidth = width;
    m_screenHeight = height;
    GLFWwindow* window = glfwCreateWindow(m_screenWidth, m_screenHeight, title.c_str(), NULL, NULL);
    if (!window) {
        glfwTerminate();
        throw std::runtime_error("Failed to create GLFW window");
    }
    glfwMakeContextCurrent(window);
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
        throw std::runtime_error("Failed to initialize GLAD");
    }
    glEnable(GL_PROGRAM_POINT_SIZE);

    m_window = window;
    glfwSwapInterval(1);
}

bool Renderer::IsRendering() const {
    return m_window && !glfwWindowShouldClose(m_window);
}

void Renderer::LoadScreen() const {
    glBindVertexArray(0);
    glfwSwapBuffers(m_window);
    glfwPollEvents();
}

void Renderer::ClearScreen() const {
    glClear(GL_DEPTH_BUFFER_BIT);
    auto [r, g, b, a] = m_backgroundColor.toScalar();
    glClearColor(r, g, b, a);
    glClear(GL_COLOR_BUFFER_BIT);
}

void Renderer::SetBackground(const Color& color) {
    m_backgroundColor = color;
}

std::string Renderer::ReadShaderSource(const char* filePath) {
    std::string content;
    std::ifstream fileStream(filePath, std::ios::in);
    std::string line;
    while (std::getline(fileStream, line)) {
        content.append(line + "\n");
    }
    fileStream.close();
    return content;
}

void Renderer::PrintShaderLog(GLuint shader) {
    Integer logLength, buffer;
    char* log;
    glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &logLength);
    if (logLength > 0) {
        log = new char[logLength];
        glGetShaderInfoLog(shader, logLength, &buffer, log);
        std::cout << "Shader Info Log: " << log << std::endl;
        delete[] log;
    }
}

void Renderer::PrintProgramLog(GLuint program) {
    Integer logLength, buffer;
    char* log;
    glGetProgramiv(program, GL_INFO_LOG_LENGTH, &logLength);
    if (logLength > 0) {
        log = new char[logLength];
        glGetProgramInfoLog(program, logLength, &buffer, log);
        std::cout << "Program Info Log: " << log << std::endl;
        delete[] log;
    }
}

bool Renderer::CheckOpenGLError() {
    bool hasError = false;
    GLenum glError = glGetError();
    while (glError != GL_NO_ERROR) {
        std::cout << "glError: " << glError << std::endl;
        hasError = true;
        glError = glGetError();
    }
    return hasError;
}

GLuint Renderer::CreateShader(const char* filePath, GLenum shaderType) {
    GLint status;
    GLuint shader = glCreateShader(shaderType);
    std::string source = ReadShaderSource(filePath);
    const char* c_source = source.c_str();
    glShaderSource(shader, 1, &c_source, NULL);
    glCompileShader(shader);
    glGetShaderiv(shader, GL_COMPILE_STATUS, &status);
    if (status == GL_FALSE) {
        PrintShaderLog(shader);
        throw std::runtime_error("Shader compilation failed");
    }
    return shader;
}

GLuint Renderer::CreateShaderProgram(const char* vp, const char* fp) {
    GLuint vShader = CreateShader(vp, GL_VERTEX_SHADER);
    GLuint fShader = CreateShader(fp, GL_FRAGMENT_SHADER);

    GLint status;
    GLuint program = glCreateProgram();
    glAttachShader(program, vShader);
    glAttachShader(program, fShader);
    glLinkProgram(program);
    glGetProgramiv(program, GL_LINK_STATUS, &status);
    if (status == GL_FALSE) {
        PrintProgramLog(program);
        throw std::runtime_error("Program linking failed");
    }
    return program;
}

template <int DIM>
void Renderer::RenderPoints(const std::vector<VectorX<DIM>>& points, const VectorX<DIM>& minCoord, const VectorX<DIM>& maxCoord, const Scalar& pointSize, const Color& color) const {
    if (points.empty()) return;

    GLuint program = 0;
    if constexpr (DIM == 2) program = CreateShaderProgram("Shader/Scene2D_vertex.glsl", "Shader/Scene2D_fragment.glsl");
    else if constexpr (DIM == 3) program = CreateShaderProgram("Shader/Scene2D_vertex.glsl", "Shader/Scene2D_fragment.glsl");
    else static_assert(DIM == 2 || DIM == 3, "Invalid template DIM value found, expected [2, 3].");
    glUseProgram(program);

    GLuint VAO, VBO;
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);

    // Bind and upload the data
    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, points.size() * sizeof(VectorX<DIM>), points.data(), GL_STATIC_DRAW);

    // Set vertex attribute
    glVertexAttribPointer(0, DIM, GL_FLOAT, GL_FALSE, sizeof(VectorX<DIM>), BUFFER_OFFSET(0));
    glEnableVertexAttribArray(0);

    auto [r, g, b, a] = color.toScalar();
    GLint minCoordLoc = glGetUniformLocation(program, "minCoord");
    GLint maxCoordLoc = glGetUniformLocation(program, "maxCoord");
    GLint pointSizeLoc = glGetUniformLocation(program, "pointSize");
    GLint pixelColorLoc = glGetUniformLocation(program, "pixelColor");
    if constexpr (DIM == 2) {
        glUniform2fv(minCoordLoc, 1, minCoord.data());
        glUniform2fv(maxCoordLoc, 1, maxCoord.data());
    }
    else if constexpr (DIM == 3) {
        glUniform3fv(minCoordLoc, 1, minCoord.data());
        glUniform3fv(maxCoordLoc, 1, maxCoord.data());
    }
    glUniform1f(pointSizeLoc, pointSize);
    glUniform4f(pixelColorLoc, r, g, b, a);

    // Draw the points as GL_POINTS (each vertex is a circle center; you'll need a fragment shader if you want them to appear as circles).
    glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(points.size()));

    // Unbind and clean up the buffers
    glBindVertexArray(0);
    glDeleteBuffers(1, &VBO);
    glDeleteVertexArrays(1, &VAO);
    glDeleteProgram(program);
}
template void Renderer::RenderPoints<2>(const std::vector<VectorX<2>>&, const VectorX<2>&, const VectorX<2>&, const Scalar&, const Color&) const;
template void Renderer::RenderPoints<3>(const std::vector<VectorX<3>>&, const VectorX<3>&, const VectorX<3>&, const Scalar&, const Color&) const;

template <int DIM>
void Renderer::RenderLines(const std::vector<VectorX<DIM>>& points, const std::vector<Integer>& indices, const VectorX<DIM>& minCoord, const VectorX<DIM>& maxCoord, const Scalar& lineWidth, const Color& color) const {
    if (points.empty() || indices.empty()) return;
    
    GLuint program = 0;
    if constexpr (DIM == 2) program = CreateShaderProgram("Shader/Scene2D_vertex.glsl", "Shader/Scene2D_fragment.glsl");
    else if constexpr (DIM == 3) program = CreateShaderProgram("Shader/Scene3D_vertex.glsl", "Shader/Scene3D_fragment.glsl");
    else static_assert(DIM == 2 || DIM == 3, "Invalid template DIM value found, expected [2, 3].");

    glUseProgram(program);
    glLineWidth(lineWidth);

    GLuint VAO, VBO, EBO;
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &EBO);

    // Bind and upload the data
    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, points.size() * sizeof(VectorX<DIM>), points.data(), GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(Integer), indices.data(), GL_STATIC_DRAW);

    // Set vertex attribute for 2D positions (2 floats)
    glVertexAttribPointer(0, DIM, GL_FLOAT, GL_FALSE, sizeof(VectorX<DIM>), BUFFER_OFFSET(0));
    glEnableVertexAttribArray(0);

    auto [r, g, b, a] = color.toScalar();
    GLint minCoordLoc = glGetUniformLocation(program, "minCoord");
    GLint maxCoordLoc = glGetUniformLocation(program, "maxCoord");
    GLint pointSizeLoc = glGetUniformLocation(program, "pointSize");
    GLint pixelColorLoc = glGetUniformLocation(program, "pixelColor");
    if constexpr (DIM == 2) {
        glUniform2fv(minCoordLoc, 1, minCoord.data());
        glUniform2fv(maxCoordLoc, 1, maxCoord.data());
    }
    else if constexpr (DIM == 3) {
        glUniform3fv(minCoordLoc, 1, minCoord.data());
        glUniform3fv(maxCoordLoc, 1, maxCoord.data());
    }
    glUniform1f(pointSizeLoc, lineWidth);
    glUniform4f(pixelColorLoc, r, g, b, a);

    glDrawElements(GL_LINES, indices.size(), GL_UNSIGNED_INT, BUFFER_OFFSET(0));

    // Unbind and clean up the buffers
    glBindVertexArray(0);
    glDeleteBuffers(1, &VBO);
    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &EBO);
    glDeleteProgram(program);
}
template void Renderer::RenderLines<2>(const std::vector<VectorX<2>>&, const std::vector<Integer>&, const VectorX<2>&, const VectorX<2>&, const Scalar&, const Color&) const;
template void Renderer::RenderLines<3>(const std::vector<VectorX<3>>&, const std::vector<Integer>&, const VectorX<3>&, const VectorX<3>&, const Scalar&, const Color&) const;

template <int DIM>
void Renderer::RenderParticles(const std::vector<Particle<DIM>>& particles, const VectorX<DIM>& minCoord, const VectorX<DIM>& maxCoord, const Scalar& pointSize, const Color& color) const {
    if (particles.empty()) return;

    GLuint program = 0;
    if constexpr (DIM == 2) program = CreateShaderProgram("Shader/Scene2D_vertex.glsl", "Shader/Scene2D_fragment.glsl");
    else if constexpr (DIM == 3) program = CreateShaderProgram("Shader/Scene3D_vertex.glsl", "Shader/Scene3D_fragment.glsl");
    else static_assert(DIM == 2 || DIM == 3, "Invalid template DIM value found, expected [2, 3].");
    glUseProgram(program);

    GLuint VAO, VBO;
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);

    // Bind and upload the data
    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, particles.size() * sizeof(Particle<DIM>), particles.data(), GL_STATIC_DRAW);

#ifdef _DEBUG
    std::cout << "offset is " << offsetof(Particle<DIM>, position) << std::endl;
#endif

    glVertexAttribPointer(0, DIM, GL_FLOAT, GL_FALSE, sizeof(Particle<DIM>), BUFFER_OFFSET(offsetof(Particle<DIM>, position)));
    glEnableVertexAttribArray(0);

    auto [r, g, b, a] = color.toScalar();
    GLint minCoordLoc = glGetUniformLocation(program, "minCoord");
    GLint maxCoordLoc = glGetUniformLocation(program, "maxCoord");
    GLint pointSizeLoc = glGetUniformLocation(program, "pointSize");
    GLint pixelColorLoc = glGetUniformLocation(program, "pixelColor");
    if constexpr (DIM == 2) {
        glUniform2fv(minCoordLoc, 1, minCoord.data());
        glUniform2fv(maxCoordLoc, 1, maxCoord.data());
    }
    else if constexpr (DIM == 3) {
        glUniform3fv(minCoordLoc, 1, minCoord.data());
        glUniform3fv(maxCoordLoc, 1, maxCoord.data());
    }
    glUniform1f(pointSizeLoc, pointSize);
    glUniform4f(pixelColorLoc, r, g, b, a);

    // Draw the points as GL_POINTS (each vertex is a circle center; you'll need a fragment shader if you want them to appear as circles).
    glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(particles.size()));

    // Unbind and clean up the buffers
    glBindVertexArray(0);
    glDeleteBuffers(1, &VBO);
    glDeleteVertexArrays(1, &VAO);
    glDeleteProgram(program);
}
template void Renderer::RenderParticles<2>(const std::vector<Particle<2>>&, const VectorX<2>&, const VectorX<2>&, const Scalar&, const Color&) const;
template void Renderer::RenderParticles<3>(const std::vector<Particle<3>>&, const VectorX<3>&, const VectorX<3>&, const Scalar&, const Color&) const;
