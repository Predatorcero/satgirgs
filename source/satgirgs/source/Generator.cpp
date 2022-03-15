#include <fstream>
#include <iostream>
#include <iomanip>
#include <random>
#include <functional>
#include <mutex>
#include <ios>

#include <omp.h>

#include <satgirgs/Generator.h>
#include <satgirgs/SpatialTree.h>


namespace satgirgs {

std::vector<double> generateWeights(int n, double ple, int weightSeed, bool parallel) {
    const auto threads = parallel ? std::max(1, std::min(omp_get_max_threads(), n / 10000)) : 1;
    auto result = std::vector<double>(n);

    #pragma omp parallel num_threads(threads)
    {
        const auto tid = omp_get_thread_num();
        auto gen = default_random_engine{weightSeed >= 0 ? (weightSeed+tid) : std::random_device()()};
        auto dist = std::uniform_real_distribution<>{};

        #pragma omp for schedule(static)
        for (int i = 0; i < n; ++i) {
            result[i] = std::pow((std::pow(0.5*n, -ple + 1) - 1) * dist(gen) + 1, 1 / (-ple + 1));
        }
    }

    return result;
}

std::vector<std::vector<double>> generatePositions(int n, int dimension, int positionSeed, bool parallel) {
    const auto threads = parallel ? std::max(1, std::min(omp_get_max_threads(), n / 10000)) : 1;
    auto result = std::vector<std::vector<double>>(n, std::vector<double>(dimension));

    #pragma omp parallel num_threads(threads)
    {
        const auto tid = omp_get_thread_num();
        auto gen = default_random_engine{positionSeed >= 0 ? (positionSeed+tid) : std::random_device()()};
        auto dist = std::uniform_real_distribution<>{};

        #pragma omp for schedule(static)
        for(int i=0; i<n; ++i)
            for (int d=0; d<dimension; ++d)
                result[i][d] = dist(gen);
    }

    return result;
}

std::vector<Node2D> convertToNodes(std::vector<std::vector<double>> positions, std::vector<double> weights, int indiceOffset){
    assert(positions.size() == weights.size());
    std::vector<Node2D> result;
    for(int i = 0; i < positions.size(); i++){
        result.push_back(Node2D(positions[i], weights[i], indiceOffset + i));
    }
    return result;
}

std::vector<std::pair<int, int>> generateEdges(const std::vector<Node2D> &c_nodes,
        const std::vector<Node2D> &nc_nodes, bool debugMode) {

    using edge_vector = std::vector<std::pair<int, int>>;
    edge_vector result;

    std::vector<std::pair<
            edge_vector,
            uint64_t[31] /* avoid false sharing */
    > > local_edges(omp_get_max_threads());

    constexpr auto block_size = size_t{1} << 20;

    std::mutex m;
    auto flush = [&] (const edge_vector& local) {
        std::lock_guard<std::mutex> lock(m);
        result.insert(result.end(), local.cbegin(), local.cend());
    };

    auto addEdge = [&](int u, int v, int tid) {
        auto& local = local_edges[tid].first;
        local.emplace_back(u,v);
        if (local.size() == block_size) {
            flush(local);
            local.clear();
        }
    };

    const auto num_threads = omp_get_max_threads();

    // TODO convert c_positions/nc_positions to vector of Node and remove template from Node
    // TODO does parallel work here?
    #pragma omp parallel for num_threads(num_threads)
    for(int clauseIndex = 0; clauseIndex < c_nodes.size(); clauseIndex++){
        auto cp = c_nodes[clauseIndex];
        const auto threadId = omp_get_thread_num();

        auto nearest = std::min_element(nc_nodes.begin(), nc_nodes.end(), [&](const Node2D& a, const Node2D& b) {return a.weightedDistance(cp) < b.weightedDistance(cp);});
        auto nearestIndex = std::distance(nc_nodes.begin(), nearest);
        auto secondNearest = std::min_element(nc_nodes.begin(), nc_nodes.end(), [&](const Node2D& a, const Node2D& b) {return (a != *nearest || b == *nearest) && a.weightedDistance(cp) < b.weightedDistance(cp);}); // ignore first minimum
        auto secondNearestIndex = std::distance(nc_nodes.begin(), secondNearest);

        if(debugMode){
            // add clause - non-clause edges
            //pad clauseIndex to distinguish from non-clause indices
            addEdge(nearestIndex, nc_nodes.size() + clauseIndex, threadId);
            addEdge(secondNearestIndex, nc_nodes.size() + clauseIndex, threadId); 
        } else {
            // add non-clause - non-clause edges
            addEdge(nearestIndex, secondNearestIndex, threadId);
        }
    }

    for(const auto& v : local_edges)
        flush(v.first);

    return result;
}


void saveDot(const std::vector<double> &weights, const std::vector<std::vector<double>> &positions,
             const std::vector<std::pair<int, int>> &graph, const std::string &file) {

    // TODO adapt saveDot code for new model
    std::ofstream f{file};
    if(!f.is_open())
        throw std::runtime_error{"Error: failed to open file \"" + file + '\"'};
    f << "graph girg {\n\toverlap=scale;\n\n";
    f << std::fixed;
    for (int i = 0; i < weights.size(); ++i) {
        f << '\t' << i << " [label=\""
          << std::setprecision(2) << weights[i] << std::setprecision(6)
          << "\", pos=\"";
        for (auto d = 0u; d < positions[i].size(); ++d)
            f << (d == 0 ? "" : ",") << positions[i][d];
        f << "\"];\n";
    }
    f << '\n';
    for (auto &edge : graph)
        f << '\t' << edge.first << "\t-- " << edge.second << ";\n";
    f << "}\n";
}

} // namespace satgirgs
