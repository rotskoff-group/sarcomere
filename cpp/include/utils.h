#ifndef UTILS_H
#define UTILS_H

#include <cmath>
#include <vector>
#include <map>
#include <algorithm>
#include <numeric>
#include <utility>
#include <omp.h>

namespace utils {

//------------------------------------------------------------------------------
// The vec structure: represents a 2D vector/point and provides common operations.
//------------------------------------------------------------------------------
struct vec {
    double x;
    double y;

    // Addition operator.
    vec operator+(const vec& p) const {
        return vec { x + p.x, y + p.y };
    }
    vec& operator+=(const vec& other) {
        x += other.x;
        y += other.y;
        return *this;
    }

    // Subtraction operator.
    vec operator-(const vec& p) const {
        return vec { x - p.x, y - p.y };
    }
    vec& operator-=(const vec& other) {
        x -= other.x;
        y -= other.y;
        return *this;
    }

    // Multiplication by scalar.
    vec operator*(double c) const {
        return vec { x * c, y * c };
    }
    friend vec operator*(double c, const vec& v) {
        return v * c;
    }

    // Division by scalar.
    vec operator/(double c) const {
        return vec { x / c, y / c };
    }

    // Euclidean norm.
    double norm() const {
        return std::sqrt(x*x + y*y);
    }

    // Unary minus.
    vec operator-() const {
        return vec { -x, -y };
    }

    // Wrap the coordinates according to periodic boundary conditions.
    // The box vector is assumed to contain the periodic lengths in x and y.
    void pbc_wrap(const std::vector<double>& box) {
        x = x - box[0] * std::round(x / box[0]);
        y = y - box[1] * std::round(y / box[1]);
    }

    // Euclidean distance (no periodic boundaries).
    double distance(const vec& p) const {
        return std::sqrt((x - p.x) * (x - p.x) + (y - p.y) * (y - p.y));
    }

    // Distance squared with periodic boundary conditions.
    double distance_squared(const vec& p, const std::vector<double>& box) const {
        double dx = x - p.x;
        double dy = y - p.y;
        dx = dx - box[0] * std::round(dx / box[0]);
        dy = dy - box[1] * std::round(dy / box[1]);
        return dx * dx + dy * dy;
    }
    
    // Distance computed with periodic boundary conditions.
    double distance(const vec& p, const std::vector<double>& box) const {
        return std::sqrt(distance_squared(p, box));
    }

    // Dot product.
    double dot(const vec& p) const {
        return x * p.x + y * p.y;
    }
    


    // Squared norm.
    double norm_squared() const {
        return x * x + y * y;
    }
};

//------------------------------------------------------------------------------
// Free function declarations (definitions are provided in utils.cpp)
//------------------------------------------------------------------------------

// Wrap a coordinate value using periodic boundary conditions.
double pbc_wrap(double x, double& box);

// Wrap an angle into the interval (–π, π] (or any 2π interval).
void angle_wrap(double& theta);

// Compute the dot product of two 2D vectors (given as four doubles).
double dotProduct(double x1, double y1, double x2, double y2);

bool compare_indices(const std::vector<int>& a, const std::vector<int>& b);


//------------------------------------------------------------------------------
// MoleculeConnection class declaration
//------------------------------------------------------------------------------
class MoleculeConnection {
public:
    MoleculeConnection();
    MoleculeConnection(int numA);

    // Add a connection from molecule A (index aIndex) to molecule B (index bIndex).
    void addConnection(int aIndex, int bIndex);
    // Delete a connection from molecule A to molecule B.
    void deleteConnection(int aIndex, int bIndex);
    // Delete all connections for molecule A.
    void deleteAllConnections(int aIndex);
    // Retrieve the connections for molecule A.
    const std::vector<int>& getConnections(int aIndex) const;

private:
    std::vector<std::vector<int>> connections;
};

//------------------------------------------------------------------------------
// Template functions (must be in the header)
//------------------------------------------------------------------------------

// Returns a vector of indices that sort the input vector in descending order.
template <typename T>
std::vector<size_t> sort_indices(const std::vector<T>& vec) {
    std::vector<size_t> indices(vec.size());
    for (size_t i = 0; i < indices.size(); ++i) {
        indices[i] = i;
    }
    std::sort(indices.begin(), indices.end(),
              [&vec](size_t a, size_t b) { return vec[a] > vec[b]; });
    return indices;
}

// Reduces a temporary array (with one sub-array per thread) into a target array.
template <typename T>
void reduce_array(std::vector<std::vector<T>>& temp_array, std::vector<T>& target_array) {
    #pragma omp for
    for (size_t i = 0; i < target_array.size(); ++i) {
        for (int t = 0; t < omp_get_num_threads(); ++t) {
            target_array[i] += temp_array[t][i];
        }
    }
}

} // namespace utils

#endif // UTILS_H
