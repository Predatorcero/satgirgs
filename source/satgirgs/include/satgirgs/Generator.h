#pragma once

#include <vector>
#include <string>

#include <satgirgs/satgirgs_api.h>
#include <satgirgs/Node.h>


namespace satgirgs {


/**
 * @brief
 *  The weights are sampled according to a power law distribution between [1, n)
 *
 * @param n
 *  The size of the graph. Should match with size of positions.
 * @param ple
 *  The power law exponent to sample the new weights. Should be 2.0 to h3.0.
 * @param weightSeed
 *  A seed for weight sampling. Should not be equal to the position seed.
 *
 * @return
 *  The weights according to the desired distribution.
 */
SATGIRGS_API std::vector<double> generateWeights(int n, double ple, int weightSeed, bool parallel = true);

/**
 * @brief
 *  Samples 2 dimensional coordinates for n points on a torus \f$[0,1)^d\f$.
 *
 * @param n
 *  Size of the graph.
 * @param positionSeed
 *  Seed to sample the positions.
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
 *  Offset for incides to distinguish clause and non-clause nodes
 *
 * @return
 *  Node vector with the positions and weights given.
 */
SATGIRGS_API std::vector<Node2D> convertToNodes(std::vector<std::vector<double>> positions, std::vector<double> weights, int indiceOffset = 0);

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
 *  An edge between node u and v is formed with probability \f$ \left(\frac{w_u w_v / W}{|| x_u - x_v ||^d}\right)^\alpha \f$ or 1.0 if the term exceeds 1.0.
 *
 * @param c_nodes
 *  Clause nodes.
 * @param nc_nodes
 *  Non-clause nodes.
 * @param debugMode
 *  In debug mode, output edges between clause and two closest non-clauses instead of the edges between these non-clauses.
 *
 * @return
 *  An edge list with zero based indices.
 */
SATGIRGS_API std::vector<std::pair<int,int>> generateEdges(const std::vector<Node2D>& c_nodes,
        const std::vector<Node2D> &nc_nodes, bool debugMode = false);


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
