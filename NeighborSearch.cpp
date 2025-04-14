#include "NeighborSearch.h"

template class HashNeighborSearch<2>;
template class HashNeighborSearch<3>;

#pragma Hash Neighbor Search2D

template <>
Integer HashNeighborSearch<2>::Hash(const VectorXi<2>& coord) const {
	Integer hashValue = ((p1 * coord.x()) ^ (p2 * coord.y()) ^ p3) % m_size;
	return (hashValue + m_size) % m_size;
}

#pragma endregion

#pragma Hash Neighbor Search3D

template <>
Integer HashNeighborSearch<3>::Hash(const VectorXi<3>& coord) const {
	Integer hashValue = ((p1 * coord.x()) ^ (p2 * coord.y()) ^ (p3 * coord.z())) % m_size;
	return (hashValue + m_size) % m_size;
}

#pragma endregion
