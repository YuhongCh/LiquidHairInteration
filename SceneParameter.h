#pragma once

#include "Definition.h"
#include "PressureSolver.h"
#include <nlohmann/json.hpp>

using json = nlohmann::json;

template <int DIM>
struct SolverParameter {
	Integer maxIterations;
	Scalar tolerance;
	CGSolver<DIM>::PrecondType preconditionerType;
	Integer numLevel;
	Integer numPreSmooth;
	Integer numPostSmooth;
	Integer numFinalSmooth;
	Scalar smoothingFactor;

	SolverParameter(const std::string& filename) {
		json j = json::parse(std::ifstream(filename));
		auto solver = j.at("SolverParameter");
		maxIterations = solver.at("maxIterations").get<Integer>();
		tolerance = solver.at("tolerance").get<Scalar>();
		numLevel = solver.at("numLevel").get<Integer>();
		numPreSmooth = solver.at("numPreSmooth").get<Integer>();
		numPostSmooth = solver.at("numPostSmooth").get<Integer>();
		numFinalSmooth = solver.at("numFinalSmooth").get<Integer>();
		smoothingFactor = solver.at("smoothingFactor").get<Scalar>();

		std::string type = solver.at("preconditionerType").get<std::string>();
		if (type == "None") {
			preconditionerType = CGSolver<DIM>::PrecondType::None;
		}
		else if (type == "Jacobi") {
			preconditionerType = CGSolver<DIM>::PrecondType::Jacobi;
		}
		else if (type == "Multigrid") {
			preconditionerType = CGSolver<DIM>::PrecondType::Multigrid;
		}
		else {
			throw std::runtime_error("Unknown preconditioner type");
		}
	}
};

template <int DIM>
struct LiquidParameter {
	Scalar density;
	Scalar viscosity;
	Integer numParticles;
	Scalar surfaceTension;
	Integer correctCycle;

	LiquidParameter(const std::string& filename) {
		json j = json::parse(std::ifstream(filename));
		auto liquid = j.at("LiquidParameter");
		density = liquid.at("density").get<Scalar>();
		viscosity = liquid.at("viscosity").get<Scalar>();
		numParticles = liquid.at("numParticles").get<Integer>();
		surfaceTension = liquid.at("surfaceTension").get<Scalar>();
		correctCycle = liquid.at("correctCycle").get<Integer>();
	}

};

template <int DIM>
struct HairParameter {

	struct HairInfo {
		VectorX<DIM> position;
		VectorX<DIM> velocity;
		bool isFixed;
	};;

	Scalar radius;
	Scalar youngsModulus;
	Scalar shearModulus;
	Scalar density;
	Scalar refRotation;
	Scalar viscosity;
	bool accumViscosity;

	Integer numHairs;
	Integer numParticles;
	std::vector<std::vector<HairInfo>> hairInfo;

	HairParameter(const std::string& filename) {
		json j = json::parse(std::ifstream(filename));
		auto h = j.at("HairParameter");

		radius = h.at("radius").get<Scalar>();
		youngsModulus = h.at("youngsModulus").get<Scalar>();
		shearModulus = h.at("shearModulus").get<Scalar>();
		viscosity = h.at("viscosity").get<Scalar>();       // 解析 viscosity
		density = h.at("density").get<Scalar>();
		refRotation = h.at("refRotation").get<Scalar>();
		accumViscosity = h.at("accumViscosity").get<bool>();

		// 数量信息
		numHairs = h.at("numHair").get<Integer>();
		numParticles = h.at("numParticles").get<Integer>();

		// info 数组
		const auto& infos = h.at("info");
		hairInfo.clear();
		hairInfo.reserve(infos.size());
		for (const auto& hairArray : infos) {
			std::vector<HairInfo> oneHair;
			oneHair.reserve(hairArray.size());
			for (const auto& elem : hairArray) {
				HairInfo info;
				// 读取 position、velocity
				for (int d = 0; d < DIM; ++d) {
					info.position[d] = elem.at("position")[d].get<Scalar>();
					info.velocity[d] = elem.at("velocity")[d].get<Scalar>();
				}
				info.isFixed = elem.at("isFixed").get<bool>();
				oneHair.push_back(info);
			}
			hairInfo.push_back(std::move(oneHair));
		}
	}

	// Clear hair info data to release the memory
	inline void ClearHairInfo() { hairInfo = {}; }
};



template <int DIM>
struct SceneParameter {
	Scalar dt;
	Scalar gravity;
	VectorX<DIM> minCoord;
	VectorX<DIM> maxCoord;
	VectorXi<DIM> gridDimension;

	SceneParameter(const std::string& filename) {
		json j = json::parse(std::ifstream(filename));
		dt = j.at("dt").get<Scalar>();
		gravity = j.at("gravity").get<Scalar>();
		auto gridDimension = j.at("gridDimension").get<std::vector<Integer>>();
		auto minCoord = j.at("minCoord").get<std::vector<Scalar>>();
		auto maxCoord = j.at("maxCoord").get<std::vector<Scalar>>();
		for (int i = 0; i < DIM; ++i) {
			this->gridDimension[i] = gridDimension[i];
			this->minCoord[i] = minCoord[i];
			this->maxCoord[i] = maxCoord[i];
		}
	}
};

