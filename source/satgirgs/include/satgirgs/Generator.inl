#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <random>
#include <functional>
#include <mutex>
#include <ios>
#include <tuple>
#include <set>

#include <omp.h>


namespace satgirgs {
template<unsigned int d>
std::vector<Node<d>> convertToNodes(std::vector<std::vector<double>> positions, std::vector<double> weights, int indexOffset){
    assert(positions.size() == weights.size());
    std::vector<Node<d>> result;
    for(int i = 0; i < positions.size(); i++){
        result.emplace_back(positions[i], weights[i], indexOffset + i);
    }
    return result;
}

// we could make this more efficient using a k-d-tree (refer to https://rosettacode.org/wiki/K-d_tree#C.2B.2B for this)
template<unsigned int d> std::vector<std::pair<int, int>> generateEdges(const std::vector<Node<d>> &c_nodes,
        const std::vector<Node<d>> &nc_nodes, int k, float t, int edgeSeed, bool debugMode) {

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
        if(u > v) std::swap(u, v);
        local.emplace_back(u,v);
        if (local.size() == block_size) {
            flush(local);
            local.clear();
        }
    };

    const auto num_threads = omp_get_max_threads();

    //const auto tid = omp_get_thread_num();
    const auto tid = 0;
    auto gen = std::default_random_engine{edgeSeed >= 0 ? (edgeSeed+tid) : std::random_device()()};

    //#pragma omp parallel for schedule(static), num_threads(num_threads)
    for(int clauseIndex = 0; clauseIndex < c_nodes.size(); clauseIndex++){
        auto cp = c_nodes[clauseIndex];
        const auto threadId = omp_get_thread_num();

        std::vector<Node<d>> clauseNodes;

        if (t == 0) {
            auto sortedNodes = std::vector<Node<d>>(nc_nodes.begin(), nc_nodes.end());
            std::partial_sort(sortedNodes.begin(), sortedNodes.begin() + k, sortedNodes.end(), [&cp](const Node<d>& a, const Node<d>& b) {return a.weightedDistance(cp) < b.weightedDistance(cp);});
            clauseNodes.reserve(k);
            std::copy(sortedNodes.begin(), sortedNodes.begin() + k, clauseNodes.begin());
        } else {
            std::vector<float> nodeWeights(nc_nodes.size());
            for (int i = 0; i < nc_nodes.size(); ++i) {
                const auto nodeWeight = nc_nodes[i].weight;
                nodeWeights[i] = std::pow(nodeWeight / std::pow(nc_nodes[i].distance(cp), d), 1 / t);
            }
            /*
            auto dist = std::discrete_distribution(nodeWeights.begin(), nodeWeights.end());
            std::set<int> clauseNodesSet;
            while (clauseNodesSet.size() < k) {
                clauseNodesSet.insert(dist(gen));
            }
            for (auto clauseNode : clauseNodesSet) {
                clauseNodes.push_back(nc_nodes[clauseNode]);
            }*/
            for (int i = 0; i < k; ++i) {
                auto dist = std::discrete_distribution(nodeWeights.begin(), nodeWeights.end());
                auto chosenNode = dist(gen);
                clauseNodes.push_back(nc_nodes[chosenNode]);
                nodeWeights[chosenNode] = 0;
            }
        }


        if(debugMode){
            // add non-clause - clause edges
            // offset clauseIndex by number of non-clause nodes to distinguish from non-clause indices
            for (int i = 0; i < k; ++i) {
                addEdge(clauseNodes[i].index, nc_nodes.size() + clauseIndex, threadId);
            }

        } else {
            // add non-clause - non-clause edges
            for (int i = 0; i < k; ++i) {
                for (int j = i + 1; j < k; ++j) {
                    addEdge(clauseNodes[i].index, clauseNodes[j].index, threadId);

                }
            }
        }
    }

    for(const auto& v : local_edges)
        flush(v.first);

    return result;
}


} // namespace satgirgs

#pragma clang diagnostic pop