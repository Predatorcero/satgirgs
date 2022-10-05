#pragma once

#include <vector>
#include <string>

#include <satgirgs/satgirgs_api.h>
#include <satgirgs/Node.h>


namespace satgirgs {


/**
 * @brief
 *  The weights are all equal.
 *
 * @param n
 *  The size of the graph. Should match with size of positions.
 * @return
 *  The weights according to the desired distribution.
 */
SATGIRGS_API std::vector<double> generateWeights(int n, double ple, int weightSeed, bool parallel = false);

/**
 * @brief
 *  Samples 2 dimensional coordinates for n points on a torus \f$[0,1)^d\f$.
 *
 * @param n
 *  Size of the graph.
 * @param positionSeed
 *  Seed to sample the positions. Should not be equal to the weight seed.
 *
 * @return
 *  The positions on a torus. All inner vectors have the same length.
 */
SATGIRGS_API std::vector<std::vector<double>> generatePositions(int n, int dimension, int positionSeed, bool parallel = true);

/**
 * @brief
 *  Zips positions and weights into one vector of 2-dimensional nodes. Both vectors should have the same length.
 *
 * @param positions
 *  Node positions
 * @param weights
 *  Node weights (power-law distributed)
 * @param indiceOffset
 *  Offset for incides to distinguish clause and non-clause nodes.
 *  E.g. when converting clause nodes, choose as offset the number of non-clause nodes.
 *
 * @return
 *  Node vector with the positions and weights given.
 */
template<unsigned int d=2> SATGIRGS_API std::vector<Node<d>> convertToNodes(std::vector<std::vector<double>> positions, std::vector<double> weights, int indexOffset = 0);

/**
 * @brief
 *  Removes duplicate edges and replaces them by a weighted edge (weight = number of occurences).
 *
 * @param edges
 *  Edges as pairs of nodes
 * @return
 *  Deduplicated edges as tuples (node, node, weight)
 */
SATGIRGS_API std::vector<std::tuple<int,int,int>> deduplicateEdges(std::vector<std::pair<int, int>> &edges);

/**
 * @brief
 *  Samples edges according to weights and positions.
 *  For each clause node c, an edge between the two non-clause nodes with the smallest weighted distance to c is formed.
 *
 * @param c_nodes
 *  Clause nodes.
 * @param nc_nodes
 *  Non-clause nodes.
 * @param k
 *  Literals per clause.
 * @param t
 *  The temperature.
 * @param debugMode
 *  In debug mode, output edges between clause and k closest non-clauses instead of the edges between these non-clauses.
 *
 * @return
 *  An edge list with zero based indices.
 */
template<unsigned int d=2> SATGIRGS_API std::vector<std::pair<int,int>> generateEdges(const std::vector<Node<d>>& c_nodes,
        const std::vector<Node<d>> &nc_nodes, int k, float t, int edgeSeed, bool debugMode = false);


/**
 * @brief
 *  Saves the graph in .dot format (graphviz).
 *  The weight is saved as a label and the coordinates as a position attribute for each Node.
 *
 * @param c_nodes
 *  Clause nodes
 * @param nc_nodes
 *  Non-clause nodes
 * @param graph
 *  An edge list with zero based indices (tuple: node, node, weight).
 * @param file
 *  The name of the output file.
 * @param debugMode
 *  In debug mode, also output clause nodes.
 */
SATGIRGS_API void saveDot(const std::vector<Node2D>& c_nodes, const std::vector<Node2D>& nc_nodes,
        const std::vector<std::tuple<int,int,int>> &graph, const std::string &file, bool debugMode = false);



} // namespace satgirgs

#include <satgirgs/Generator.inl>