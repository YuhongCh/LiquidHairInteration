#pragma once
#include "Definition.h"

template <typename T>
class Array2 {
public:

    Array2(Integer lenX = 1, Integer lenY = 1)
        : m_dimension(lenX, lenY), m_container(lenX * lenY) {
    }

    Array2(const Vector2i& dimension)
        : m_dimension(dimension), m_container(dimension[0] * dimension[1]) {
    }

    inline typename std::vector<T>::const_reference operator()(const Integer& xi, const Integer& yi) const {
        return m_container[xi * m_dimension[1] + yi];
    }
    inline typename std::vector<T>::reference operator()(const Integer& xi, const Integer& yi) {
        return m_container[xi * m_dimension[1] + yi];
    }
    inline const Vector2i& GetDimension() const { return m_dimension; }

    inline Integer GetSize() const { return static_cast<Integer>(m_container.size()); }

    inline void Resize(const Integer& xi, const Integer& yi) {
        m_dimension = Vector2i(xi, yi);
        m_container.resize(xi * yi);
    }

    inline void Resize(const Vector2i& dim) const { Resize(dim.x(), dim.y()); }

    inline void Fill(const T& val) {
        Integer N = GetSize();

        #pragma omp simd
        for (Integer i = 0; i < N; ++i) {
            m_container[i] = val;
        }
    }

    inline void CopyFrom(const Array2<T>& source) {
		Integer N = GetSize();
#pragma omp simd
		for (Integer i = 0; i < N; ++i) {
			m_container[i] = source.m_container[i];
		}
    }

    inline static T Reduce(const Array2& arr1, const Array2& arr2) {
        assert(arr1.m_dimension == arr2.m_dimension);
        Integer N = arr1.GetSize();

        T sum = T();
        #pragma omp parallel for reduction(+:sum) schedule(static)
        for (Integer i = 0; i < N; ++i) {
            sum += arr1.m_container[i] * arr2.m_container[i];
        }
        return sum;
    }

private:
    Vector2i m_dimension;
    std::vector<T> m_container;
};

template <typename T>
class Array3 {
public:
    Array3(Integer lenX = 1, Integer lenY = 1, Integer lenZ = 1)
        : m_dimension(lenX, lenY, lenZ), m_container(lenX * lenY * lenZ) {
    }

    Array3(const Vector3i& dimension)
        : m_dimension(dimension), m_container(dimension[0] * dimension[1] * dimension[2]) {
    }

    inline typename std::vector<T>::const_reference operator()(const Integer& xi, const Integer& yi, const Integer& zi) const {
        return m_container[(xi * m_dimension.y() + yi) * m_dimension.z() + zi];
    }
    inline typename std::vector<T>::reference operator()(const Integer& xi, const Integer& yi, const Integer& zi) {
        return m_container[(xi * m_dimension.y() + yi) * m_dimension.z() + zi];
    }
    inline const Vector3i& GetDimension() const { return m_dimension; }

    inline Integer GetSize() const { return static_cast<Integer>(m_container.size()); }

    inline void Resize(const Integer& xi, const Integer& yi, const Integer& zi) {
        m_dimension = Vector3i(xi, yi, zi);
        m_container.resize(xi * yi * zi);
    }

    inline void Resize(const Vector3i& dim) { Resize(dim.x(), dim.y(), dim.z()); }

    inline void Fill(const T& val) {
        Integer N = GetSize();
        #pragma omp simd
        for (Integer i = 0; i < N; ++i) {
            m_container[i] = val;
        }
    }

    inline void CopyFrom(const Array3<T>& source) {
        Integer N = GetSize();
#pragma omp simd
        for (Integer i = 0; i < N; ++i) {
            m_container[i] = source.m_container[i];
        }
    }

    inline static T Reduce(const Array3& arr1, const Array3& arr2) {
        assert(arr1.m_dimension == arr2.m_dimension);
        Integer N = arr1.GetSize();

        T sum = T();
#pragma omp parallel for reduction(+:sum) schedule(static)
        for (Integer i = 0; i < N; ++i) {
            sum += arr1.m_container[i] * arr2.m_container[i];
        }
        return sum;
    }

private:
    Vector3i m_dimension;
    std::vector<T> m_container;
};
