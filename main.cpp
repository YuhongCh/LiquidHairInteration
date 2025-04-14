#include "Definition.h"
#include "Renderer.h"
#include "Scene.h"
#include "LiquidModel.h"
#include "Boundary.h"

#ifdef _DEBUG
bool pauseSimulation = false;

void OnLeftMouseClicked(GLFWwindow* window, int button, int action, int mods) {
	if (button == GLFW_MOUSE_BUTTON_LEFT) {
		if (action == GLFW_PRESS) {
			int width, height;
			double xpos, ypos;
			glfwGetCursorPos(window, &xpos, &ypos);
			glfwGetWindowSize(window, &width, &height);

			std::cout << "Left mouse button pressed at (" << xpos / width << ", " << ypos / height << ")" << std::endl;
		}
	}
}

void OnSpaceClicked(GLFWwindow* window, int key, int scancode, int action, int mods) {
	if (key == GLFW_KEY_SPACE && action == GLFW_PRESS) {
		pauseSimulation = !pauseSimulation;
		if (pauseSimulation) {
			std::cout << "Space button pressed, Simulation Paused " << std::endl;
		}
		else {
			std::cout << "Space button pressed, Simulation Continue to play " << std::endl;
		}
	}
}

#endif

int main(int argc, char* argv[]) {

	Renderer& rend = Renderer::GetInstance();
	rend.InstantiateWindow("OpenGL Window", 480, 480);
	rend.SetBackground(Color::White());
	
#ifdef _DEBUG
	rend.SetMouseCallback(OnLeftMouseClicked);
	rend.SetKeyboardCallback(OnSpaceClicked);
#endif

	SceneParameter<DIMENSION> params = SceneParameter<DIMENSION>::Default();
	Scene<DIMENSION> scene(params);
	LiquidModel<DIMENSION>& model = scene.GetLiquidModel();

#ifdef _DEBUG
#if DIMENSION == 2
	std::vector<Vector2> points = {Vector2(-0.9, -0.9), Vector2(0.9, -0.9), Vector2(0.9, 0.9), Vector2(-0.9, 0.9)};
	std::vector<Integer> indices = { 0, 1, 1, 2, 2, 3, 3, 0 };
#else
	std::vector<Vector3> points = { Vector3(-0.9, -0.9, -0.9), Vector3(-0.9, -0.9, 0.9), Vector3(-0.9, 0.9, -0.9), Vector3(-0.9, 0.9, 0.9),
									Vector3(0.9, -0.9, -0.9), Vector3(0.9, -0.9, 0.9), Vector3(0.9, 0.9, -0.9), Vector3(0.9, 0.9, 0.9) };
	std::vector<Integer> indices = { 0, 1, 1, 5, 5, 4, 4, 0, 2, 3, 3, 7, 7, 6, 6, 2, 0, 2, 1, 3, 4, 6, 5, 7 };
#endif
#endif

	model.Init();
	Scalar dt = 0.01;


	while (rend.IsRendering()) {
		rend.ClearScreen();

#ifdef _DEBUG
		if (pauseSimulation) {
			rend.LoadScreen();
			continue;
		}
		std::cout << "Start step with dt=" << dt << std::endl;
#endif

		model.PreStep(dt);
		model.Step(dt);
		model.PostStep(dt);
		dt = model.GetCFL();

		rend.RenderParticles(model.GetParticleSystem().GetParticles(), params.minCoord, params.maxCoord, 5.0, Color::Blue());

#ifdef _DEBUG
		rend.RenderLines(points, indices, params.minCoord, params.maxCoord, 5.0, Color::Black());
#endif

		rend.LoadScreen();
	}
	
}