#include "Kappa.h"

#pragma region template declaration
template class Kappa<3>;
template class GradKappa<3>;
template class HessKappa<3>;

#pragma endregion

#pragma region Kappa

template <>
Kappa<3>::Kappa(const std::vector<VectorX<3>>& kb, const std::vector<VectorX<3>>& matFrame1, const std::vector<VectorX<3>>& matFrame2)
	: m_kb(kb), m_matFrame1(matFrame1), m_matFrame2(matFrame2), m_size(static_cast<Integer>(kb.size())) {
	m_kappa.resize(m_size);
}

template <>
void Kappa<3>::Compute() {
#pragma omp parallel for
	for (Integer index = 0; index < m_size; ++index) {
		m_kappa[index] = Vector4(
			m_matFrame2[index].dot(m_kb[index]), -m_matFrame1[index].dot(m_kb[index]),
			m_matFrame2[index + 1].dot(m_kb[index]), -m_matFrame1[index + 1].dot(m_kb[index])
		);
	}
}

template <>
void Kappa<3>::CopyFrom(const Kappa<3>& other) {
#ifdef _DEBUG
	ASSERT_MSG(m_size == other.m_size, "Size mismatch!");
#endif
	for (Integer i = 0; i < m_size; ++i) {
		m_kappa[i] = other.m_kappa[i];
	}
}

template <>
void Kappa<3>::Print() const {
	for (Integer i = 0; i < m_size; ++i) {
		std::cout << "Kappa[" << i << "]: " << m_kappa[i].transpose() << "\n";
	}
}

#pragma endregion

#pragma region GradKappa
template <>
GradKappa<3>::GradKappa(const std::vector<VectorX<3>>& tangent,
	const std::vector<Scalar>& length,
	const std::vector<VectorX<3>>& kb,
	const std::vector<VectorX<3>>& matFrame1,
	const std::vector<VectorX<3>>& matFrame2,
	const std::vector<Vector4>& kappa)
	: m_tangent(tangent), m_length(length), m_kb(kb), 
	m_matFrame1(matFrame1), m_matFrame2(matFrame2), m_kappa(kappa), m_size(static_cast<Integer>(kappa.size())) {
	m_gradKappa.resize(m_size);
}

void GradKappa<3>::Compute() {
#pragma omp parallel for
    for (Integer index = 0; index < m_size; ++index) {
		const Vector3& ts = m_tangent[index];
        const Vector3& te = m_tangent[index + 1];
		const Scalar& ls = m_length[index];
		const Scalar& le = m_length[index + 1];
        const Vector3& m1s = m_matFrame1[index];
        const Vector3& m1e = m_matFrame1[index + 1];
        const Vector3& m2s = m_matFrame2[index];
        const Vector3& m2e = m_matFrame2[index + 1];
        const Vector4& kappa = m_kappa[index];
		MatrixX<11, 4>& gradKappa = m_gradKappa[index];

        Scalar chi = MathUtils::Clamp(1.0 + ts.dot(te), EPSILON, INF);
		Scalar inv_chi = 1.0 / chi;

		Vector3 tilde_t = inv_chi * (ts + te);
		Vector3 Dks1Des = (-kappa[0] * tilde_t + 2 * te.cross(m2s) * inv_chi) / ls;
		Vector3 Dks1Dee = (-kappa[0] * tilde_t - 2 * ts.cross(m2s) * inv_chi) / le;
		Vector3 Dks2Des = (-kappa[1] * tilde_t - 2 * te.cross(m1s) * inv_chi) / ls;
		Vector3 Dks2Dee = (-kappa[1] * tilde_t + 2 * ts.cross(m1s) * inv_chi) / le;
		Vector3 Dke1Des = (-kappa[2] * tilde_t + 2 * te.cross(m2e) * inv_chi) / ls;
		Vector3 Dke1Dee = (-kappa[2] * tilde_t - 2 * ts.cross(m2e) * inv_chi) / le;
		Vector3 Dke2Des = (-kappa[3] * tilde_t - 2 * te.cross(m1e) * inv_chi) / ls;
		Vector3 Dke2Dee = (-kappa[3] * tilde_t + 2 * ts.cross(m1e) * inv_chi) / le;

		gradKappa.block<3, 1>(0, 0) = -Dks1Des;
		gradKappa.block<3, 1>(4, 0) = Dks1Des - Dks1Dee;
		gradKappa.block<3, 1>(8, 0) = Dks1Dee;
		gradKappa.block<3, 1>(0, 1) = -Dks2Des;
		gradKappa.block<3, 1>(4, 1) = Dks2Des - Dks2Dee;
		gradKappa.block<3, 1>(8, 1) = Dks2Dee;

		gradKappa.block<3, 1>(0, 2) = -Dke1Des;
		gradKappa.block<3, 1>(4, 2) = Dke1Des - Dke1Dee;
		gradKappa.block<3, 1>(8, 2) = Dke1Dee;
		gradKappa.block<3, 1>(0, 3) = -Dke2Des;
		gradKappa.block<3, 1>(4, 3) = Dke2Des - Dke2Dee;
		gradKappa.block<3, 1>(8, 3) = Dke2Dee;

		gradKappa(3, 0) = -m_kb[index].dot(m1s);
		gradKappa(7, 0) = 0.0;
		gradKappa(3, 1) = -m_kb[index].dot(m2s);
		gradKappa(7, 1) = 0.0;
		gradKappa(3, 2) = 0.0;
		gradKappa(7, 2) = -m_kb[index].dot(m1e);
		gradKappa(3, 3) = 0.0;
		gradKappa(7, 3) = -m_kb[index].dot(m2e);
    }
}

void GradKappa<3>::CopyFrom(const GradKappa<3>& other) {
#ifdef _DEBUG
	ASSERT_MSG(m_size == other.m_size, "Size mismatch!");
#endif
	for (Integer i = 0; i < m_size; ++i) {
		m_gradKappa[i] = other.m_gradKappa[i];
	}
}

void GradKappa<3>::Print() const {
	for (Integer i = 0; i < m_size; ++i) {
		std::cout << "GradKappa[" << i << "]: " << m_gradKappa[i].transpose() << "\n";
	}
}

#pragma endregion

#pragma region HessKappa
template <>
HessKappa<3>::HessKappa(const std::vector<VectorX<3>>& tangent,
						const std::vector<Scalar>& length,
						const std::vector<VectorX<3>>& kb,
						const std::vector<VectorX<3>>& matFrame1,
						const std::vector<VectorX<3>>& matFrame2,
						const std::vector<Vector4>& kappa)
	: m_tangent(tangent), m_length(length), m_kb(kb),
	m_matFrame1(matFrame1), m_matFrame2(matFrame2), m_kappa(kappa), m_size(static_cast<Integer>(kappa.size())) {
	m_hessKappa.resize(m_size);
}

void HessKappa<3>::CopyFrom(const HessKappa<3>& other) {
#ifdef _DEBUG
	ASSERT_MSG(m_size == other.m_size, "Size mismatch!");
#endif
	for (Integer i = 0; i < m_size; ++i) {
		m_hessKappa[i] = other.m_hessKappa[i];
	}
}

void HessKappa<3>::Print() const {
	for (Integer i = 0; i < m_size; ++i) {
		std::cout << "HessKappa[" << i << "]: " << m_hessKappa[i].transpose() << "\n";
	}
}

// TODO: finish this HessKappa computation
void HessKappa<3>::Compute() {
//#pragma omp parallel for
//    for (Integer i = 0; i < m_size; ++i) {
//        const Scalar& ls = m_length[i];
//        const Scalar& le = m_length[i + 1];
//        const Vector3& ts = m_tangent[i];
//        const Vector3& te = m_tangent[i + 1];
//        const Vector3& m1s = m_matFrame1[i];
//        const Vector3& m2s = m_matFrame2[i];
//        const Vector3& m1e = m_matFrame1[i + 1];
//        const Vector3& m2e = m_matFrame2[i + 1];
//        const Vector4& kappa = m_kappa[i];
//        const Vector3& kb = m_kb[i];
//		MatrixX<11, 11>& hessKappa = m_hessKappa[i];
//
//		Matrix3 I = Matrix3::Identity();
//		Scalar chi = MathUtils::Clamp(1.0 + ts.dot(te), EPSILON, INF);
//		Scalar inv_chi = 1.0 / chi;
//		Vector3 tilde_t = inv_chi * (ts + te);
//		Vector3 tilde_d1s = inv_chi * (m1s + m1s);
//		Vector3 tilde_d1e = inv_chi * (m1e + m1e);
//		Vector3 tilde_d2s = inv_chi * (m2s + m2s);
//		Vector3 tilde_d2e = inv_chi * (m2e + m2e);
//
//		Vector3 Dk0sDs = (-kappa[0] * tilde_t + te.cross(tilde_d2s)) / ls;
//		Vector3 Dk0sDe = (-kappa[0] * tilde_t - ts.cross(tilde_d2s)) / le;
//		Vector3 Dk1sDs = (-kappa[1] * tilde_t - te.cross(tilde_d1s)) / ls;
//		Vector3 Dk1sDe = (-kappa[1] * tilde_t + ts.cross(tilde_d1s)) / le;
//		Vector3 Dk0eDs = (-kappa[2] * tilde_t + te.cross(tilde_d2e)) / ls;
//		Vector3 Dk0eDe = (-kappa[2] * tilde_t - ts.cross(tilde_d2e)) / le;
//		Vector3 Dk1eDs = (-kappa[3] * tilde_t - te.cross(tilde_d1e)) / ls;
//		Vector3 Dk1eDe = (-kappa[3] * tilde_t + ts.cross(tilde_d1e)) / le;
//
//		Vector3 DchiDs = (I - ts * ts.transpose()) * te / ls;
//		Vector3 DchiDe = (I - te * te.transpose()) * ts / le;
//		Vector3 DteDe = (I - te * te.transpose()) * ts / le;
//
//		Vector3 DttDs = (I - ts * ts.transpose() - tilde_t * (I - ts * ts.transpose()).transpose()) * te / (chi * ls);
//		Vector3 DttDe = (I - te * te.transpose() - tilde_t * (I - te * te.transpose()).transpose()) * ts / (chi * le);
//
//		
//    }
}




#pragma endregion