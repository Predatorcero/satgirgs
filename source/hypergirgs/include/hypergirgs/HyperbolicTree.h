#pragma once

#include <vector>
#include <utility>
#include <random>

#include <hypergirgs/AngleHelper.h>
#include <hypergirgs/RadiusLayer.h>
#include <hypergirgs/Point.h>


namespace hypergirgs {


template <typename EdgeCallback>
class HyperbolicTree
{
public:

    HyperbolicTree(std::vector<double>& radii, std::vector<double>& angles, double T, double R, EdgeCallback& edgeCallback, bool profile = false);

    void generate(int seed);

protected:


    void visitCellPair(unsigned int cellA, unsigned int cellB, unsigned int level);

    void sampleTypeI(unsigned int cellA, unsigned int cellB, unsigned int level, unsigned int i, unsigned int j);

    void sampleTypeII(unsigned int cellA, unsigned int cellB, unsigned int level, unsigned int i, unsigned int j);

    unsigned int partitioningBaseLevel(double r1, double r2); // takes lower bound on radius for two layers


protected:
    EdgeCallback& m_edgeCallback;
    const bool m_profile;

    const size_t m_n; ///< number of nodes

    const double m_coshR; ///< = cosh(R)

    const double m_T; ///< temperature
    const double m_R; ///< radius

    unsigned int m_layers; ///< number of layers
    unsigned int m_levels; ///< number of levels

    std::shared_ptr<Point[]> m_points;                     ///< vector of points
    std::shared_ptr<unsigned int[]> m_first_point_in_cell; ///< prefix sum into m_points
    std::vector<RadiusLayer> m_radius_layers;

    std::vector<std::vector<std::pair<unsigned int, unsigned int> > > m_layer_pairs;

    hypergirgs::default_random_engine m_gen; ///< random generator
    std::uniform_real_distribution<> m_dist; ///< random distribution

#ifndef NDEBUG
    long long m_type1_checks; ///< number of node pairs per thread that are checked via a type 1 check
    long long m_type2_checks; ///< number of node pairs per thread that are checked via a type 2 check
#endif // NDEBUG
};

template <typename EdgeCallback>
inline HyperbolicTree<EdgeCallback> makeHyperbolicTree(std::vector<double>& radii, std::vector<double>& angles, double T, double R, EdgeCallback& edgeCallback, bool profile = false) {
    return {radii, angles, T, R, edgeCallback, profile};
}

} // namespace hypergirgs

#include <hypergirgs/HyperbolicTree.inl>
