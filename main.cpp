#include "Definition.h"
#include "Renderer.h"
#include "Scene.h"
#include "LiquidModel.h"
#include "HairModel.h"
#include "Boundary.h"

int main(int argc, char* argv[]) {

	Renderer& rend = Renderer::GetInstance();
	rend.InstantiateWindow("OpenGL Window", 480, 480);
	rend.SetBackground(Color::White());

	Scene<DIMENSION> scene("Config/PureLiquid3D.json");
	LiquidModel<DIMENSION>& liquidModel = scene.GetLiquidModel();
	HairModel<DIMENSION>& hairModel = scene.GetHairModel();

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

	//scene.PlaySimulation("LiquidData.csv");

	liquidModel.Init();
	Scalar dt = 0.01;

	//Scalar currTime = 0.0, endTime = 20.0;

	//liquidModel.GetParticleSystem().WriteToFile("LiquidData.csv", false);
	//while (currTime < endTime) {
	//	liquidModel.PreStep(dt);
	//	liquidModel.Step(dt);
	//	liquidModel.PostStep(dt);
	//	dt = liquidModel.GetCFL();
	//	liquidModel.GetParticleSystem().WriteToFile("LiquidData.csv", true);
	//	currTime += dt;
	//}


	while (rend.IsRendering()) {
		rend.ClearScreen();

		liquidModel.PreStep(dt);
		liquidModel.Step(dt);
		liquidModel.PostStep(dt);
		dt = liquidModel.GetCFL();

		rend.RenderParticles(liquidModel.GetParticleSystem().GetParticles(), scene.GetMinCoord(), scene.GetMaxCoord(), 5.0, Color::Blue());
		//hairModel.Render();
#ifdef _DEBUG
		rend.RenderLines(points, indices, scene.GetMinCoord(), scene.GetMaxCoord(), 5.0, Color::Black());
#endif

		rend.LoadScreen();
	}
	
}