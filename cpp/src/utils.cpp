#include "utils.h"
#include <algorithm>  // For std::find
#include <cmath>      // For std::sqrt, std::round, and M_PI
#include <stdexcept>


namespace utils {

//------------------------------------------------------------------------------
// Definitions of free functions
//------------------------------------------------------------------------------

// Parse PBC mode from string
PBCMask parse_pbc_mask(const std::string& mode) {
    if (mode == "periodic") {
        return PBCMask{true, true};
    } else if (mode == "mixed") {
        return PBCMask{true, false};
    } else {
        throw std::invalid_argument("Invalid PBC mode: " + mode + ". Use 'periodic' or 'mixed'.");
    }
}

void angle_wrap(double& theta) {
    theta = theta - 2 * M_PI * std::round(theta / (2 * M_PI));
}

bool compare_indices(const std::vector<int>& a, const std::vector<int>& b) {
    if (a.size() != b.size()) {
        return false;
    }
    std::vector<int> sorted_a = a;
    std::vector<int> sorted_b = b;
    std::sort(sorted_a.begin(), sorted_a.end());
    std::sort(sorted_b.begin(), sorted_b.end());
    return std::equal(sorted_a.begin(), sorted_a.end(), sorted_b.begin());
}

//------------------------------------------------------------------------------
// Definitions of MoleculeConnection member functions
//------------------------------------------------------------------------------

MoleculeConnection::MoleculeConnection() {
    // Default constructor.
}

MoleculeConnection::MoleculeConnection(int numA) : connections(numA) {
    // Initialize each molecule A's connection vector.
    for (int i = 0; i < numA; i++) {
        connections[i] = std::vector<int>();
    }
}

void MoleculeConnection::addConnection(int aIndex, int bIndex) {
    if (aIndex >= 0 && aIndex < static_cast<int>(connections.size())) {
        // Avoid adding duplicate connections.
        if (std::find(connections[aIndex].begin(), connections[aIndex].end(), bIndex) == connections[aIndex].end()) {
            connections[aIndex].push_back(bIndex);
        }
    }
}

void MoleculeConnection::deleteConnection(int aIndex, int bIndex) {
    if (aIndex >= 0 && aIndex < static_cast<int>(connections.size())) {
        auto& connList = connections[aIndex];
        auto it = std::find(connList.begin(), connList.end(), bIndex);
        if (it != connList.end()) {
            connList.erase(it);
        }
    }
}

void MoleculeConnection::deleteAllConnections(int aIndex) {
    if (aIndex >= 0 && aIndex < static_cast<int>(connections.size())) {
        connections[aIndex].clear();
    }
}

const std::vector<int>& MoleculeConnection::getConnections(int aIndex) const {
    if (aIndex >= 0 && aIndex < static_cast<int>(connections.size())) {
        return connections[aIndex];
    } else {
        static const std::vector<int> emptyList;
        return emptyList;
    }
}

} // namespace utils
