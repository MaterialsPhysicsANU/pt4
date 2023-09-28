#ifndef VEC_H
#define VEC_H

#include<array>
#include<vector>
#include<algorithm>
#include<stdexcept>
#include<cmath>
#include<iostream>

////////////////////////////////////////////////////////////////////////////////////////
    // array alias types
///////////////////////////////////////////////////////////////////////////////////////
typedef std::array<double,3> double3;
typedef std::array<float,3> float3;
typedef std::array<int,3> int3;
typedef std::array<size_t,3> size_t3;

////////////////////////////////////////////////////////////////////////////////////////
    // array operator
///////////////////////////////////////////////////////////////////////////////////////

template<typename T, typename S, std::size_t N>
inline std::array<T,N> array_static_cast(const std::array<S,N>& p) {
    std::array<T,N> result;

    for (std::size_t i = 0; i < N; ++i) {
        result[i] = static_cast<T>(p[i]);
    }

    return result;
}

template<typename T, std::size_t N>
inline void operator+=(std::array<T,N>& p, const std::array<T,N>& q) {
    for (std::size_t i = 0; i < N; ++i) {
        p[i] += q[i];
    }
}

template<typename T, std::size_t N>
inline std::array<T,N> operator+(std::array<T,N> p, const std::array<T,N>& q) {
    p += q; return p;
}

template<typename T, std::size_t N>
inline void operator+=(std::array<T,N>& p, const T q) {
    for (std::size_t i = 0; i < N; ++i) {
        p[i] += q;
    }
}

template<typename T, std::size_t N>
inline std::array<T,N> operator+(std::array<T,N> p, const T q) {
    p += q; return p;
}

template<typename T, std::size_t N>
inline std::array<T,N> operator+(const T q, std::array<T,N> p) {
    p += q; return p;
}

template<typename T, std::size_t N>
inline void operator-=(std::array<T,N>& p, const std::array<T,N>& q) {
    for (std::size_t i = 0; i < N; ++i) {
        p[i] -= q[i];
    }
}

template<typename T, std::size_t N>
inline std::array<T,N> operator-(std::array<T,N> p, const std::array<T,N>& q) {
    p -= q; return p;
}

template<typename T, std::size_t N>
inline void operator*=(std::array<T, N>& p, const std::array<T, N>& q) {
    for (std::size_t i = 0; i < N; ++i) {
        p[i] *= q[i];
    }
}

template<typename T, std::size_t N>
inline std::array<T, N> operator*(std::array<T, N> p, const std::array<T, N>& q) {
    p *= q; return p;
}

template<typename T, std::size_t N>
inline void operator*=(std::array<T, N> &p, const T v) {
    for (std::size_t i = 0; i < N; ++i) {
        p[i] *= v;
    }
}

template<typename T, std::size_t N>
inline std::array<T, N> operator*(std::array<T, N> q, const T v) {
    q *= v; return q;
}

template<typename T, std::size_t N>
inline std::array<T, N> operator*(const T v, std::array<T, N> q) {
    q *= v; return q;
}

template<typename T, std::size_t N>
inline void operator/=(std::array<T, N>& p, const std::array<T, N>& q) {
    for (std::size_t i = 0; i < N; ++i) {
        p[i] /= q[i];
    }
}

template<typename T, std::size_t N>
inline std::array<T, N> operator/(std::array<T, N> p, const std::array<T, N>& q) {
    p /= q; return p;
}

template<typename T, std::size_t N>
inline void operator/=(std::array<T, N> &p, const T v) {
    for (std::size_t i = 0; i < N; ++i) {
        p[i] /= v;
    }
}

template<typename T, std::size_t N>
inline std::array<T, N> operator/(std::array<T, N> q, const T v) {
    q /= v; return q;
}

template<typename T, std::size_t N>
inline void operator/=(const T v, std::array<T, N> &p) {
    for (std::size_t i = 0; i < N; ++i) {
        p[i] = v/p[i];
    }
}

template<typename T, std::size_t N>
inline std::array<T, N> operator/(const T v, std::array<T, N> q) {
    v /= q; return q;
}

template<typename T, std::size_t N>
inline T dot(std::array<T, N> p, const std::array<T, N>& q) {
    T inner_product = T();
    for(std::size_t i = 0; i < N; ++i){
        inner_product += p[i]*q[i];
    }
    return inner_product;
}

////////////////////////////////////////////////////////////////////////////////////////
    // vector operator
///////////////////////////////////////////////////////////////////////////////////////


template<typename T>
inline void operator+=(std::vector<T>& p, const std::vector<T>& q) {
	const std::size_t N = p.size();
    for (std::size_t i = 0; i < N; ++i) {
        p[i] += q[i];
    }
}

template<typename T>
inline std::vector<T> operator+(std::vector<T> p, const std::vector<T>& q) {
    p += q; return p;
}

template<typename T>
inline void operator-=(std::vector<T>& p, const std::vector<T>& q) {
	const std::size_t N = p.size();
    for (std::size_t i = 0; i < N; ++i) {
        p[i] -= q[i];
    }
}

template<typename T>
inline std::vector<T> operator-(std::vector<T> p, const std::vector<T>& q) {
    p -= q; return p;
}

template <typename T>
inline void operator*=(std::vector<T> &p, const T v) {
	const std::size_t N = p.size();
    for (std::size_t i = 0; i < N; ++i) {
        p[i] *= v;
    }
}

template <typename T, std::size_t N>
inline void operator*=(std::vector<std::array<T,N>> &p, const std::vector<T>& v) {
	const std::size_t s = p.size();
    for (std::size_t i = 0; i < s; ++i) {
        for (std::size_t d = 0; d < N; ++d) {
            p[i][d] *= v[i];
        }
    }
}

template <typename T>
inline void operator*=(std::vector<T> &p, const std::vector<T> &q) {
	const std::size_t N = p.size();
    for (std::size_t i = 0; i < N; ++i) {
        p[i] *= q[i];
    }
}

template<typename T>
inline std::vector<T> operator*(const T v, std::vector<T> q) {
    q *= v; return q;
}

template<typename T>
inline std::vector<T> operator*(std::vector<T> q, const T v) {
    q *= v; return q;
}

template <typename T>
inline void operator/=(std::vector<T>& p, const std::vector<T>& q) {
	const std::size_t N = p.size();
    for (std::size_t i = 0; i < N; ++i) {
        p[i] /= q[i];
    }
}

////////////////////////////////////////////////////////////////////////////////////////
    // other methods
///////////////////////////////////////////////////////////////////////////////////////

template <typename T, std::size_t N>
inline T max(const std::array<T,N>& p){
    return std::max_element(p.begin(), p.end())[0];
}

template <typename T, std::size_t N>
inline T min(const std::array<T,N>& p){
    return std::min_element(p.begin(), p.end())[0];
}


template <typename T>
inline void fill(std::vector<T> &p, const T v) {
    for (std::size_t i = 0; i < p.size(); i++) {
        p[i] = v;
    }
}

template<typename S, typename T, std::size_t N>
inline std::array<S, N> floor(const std::array<T, N>& p, std::array<T, N>& q) {
    std::array<S, N> result;

    for (size_t i = 0; i < N; i++) {
        result[i] = (S)floor(p[i]);
        q[i] = p[i] - result[i];
    }
    return result;
}

template<typename S, typename T, std::size_t N>
inline std::array<S, N> floor(const std::array<T, N>& p) {
    std::array<S, N> result;
    for (size_t i = 0; i < N; i++) {
        result[i] = (S)floor(p[i]);
    }
    return result;
}

template <typename T>
inline void clamp(std::vector<T> &p, const T min, const T max) {
    for (std::size_t i = 0; i < p.size(); i++) {
        p[i] = std::clamp(p[i], min, max);
    }
}

//the magnitude of a vector
template <typename T>
inline typename T::value_type mag(const T &p) {
    typename T::value_type sumsquare = 0;
    for (std::size_t i = 0; i < p.size(); ++i) {
        sumsquare += p[i]*p[i];
    }
    return std::sqrt(sumsquare);
}

//the square of the magnitude of a vector
template <typename T>
inline typename T::value_type mag_sq(const T &p) {
    typename T::value_type sumsquare = 0;
    for (std::size_t i = 0; i < p.size(); ++i) {
        sumsquare += p[i]*p[i];
    }
    return sumsquare;
}

template <typename T>
inline T abs(const T &p) {
    auto q = p;
    for (std::size_t i = 0; i < q.size(); ++i) {
        q[i] = std::abs(p[i]);
    }
    return q;
}

template <typename T>
inline T unit(const T &p) {
    auto q = p;
    auto mg = mag(p);
    for (std::size_t i = 0; i < q.size(); ++i) {
        q[i] /= mg;
    }
    return q;
}

template<typename T, std::size_t N>
inline T product(const std::array<T, N>& p, const size_t M = N) {
    T prod = T(1);

    for (size_t i = 0; i < M; i++) {
        prod *= p[i];
    }


    return prod;
}

template<typename T, std::size_t N>
inline T flatten(const std::array<T, N>& index, const std::array<T, N>& size) {
    T sum = index[0];

    for (size_t i = 1; i < N; i++) {
        T prod = size[0];
        for (size_t j = 1; j < i; j++) {
            prod *= size[j];
        }
    sum += index[i]*prod;
    }

    return sum;
}

template<typename T, std::size_t N>
inline std::array<T, N> unflatten(T index, const std::array<T, N>& size) {
    std::array<T,N> vec_index;

    for(size_t i = N-1; i < N; --i){
        T stride = T(1);
        for(size_t j = 0; j < i; ++j){stride *= size[j];}
        T axis_component = index/stride;
        index -= axis_component*stride;
        vec_index[i] = axis_component;
    }

    return vec_index;
}

template<typename T, std::size_t N>
std::ostream& operator<<(std::ostream& os, const std::array<T, N>& arr) {
    os << '{';
    for (std::size_t i = 0; i < N; ++i) {
        os << arr[i] << ", ";
    }
    os << '}';
    return os;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& p) {
    const size_t N = p.size();
    os << '{';
    for (std::size_t i = 0; i < N; ++i) {
        os << p[i] << ", ";
    }
    os << '}';
    return os;
}

#endif 
