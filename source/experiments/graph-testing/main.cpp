
#include <iostream>
#include <cassert>
#include <limits>

#include <omp.h>
#include "satgirgs/Generator.h"


using namespace std;


void measure(int n, int m, int k, float t, int threads, int seed, int plot) {

    omp_set_num_threads(threads);
    assert(threads == omp_get_max_threads());

    auto ncseed = seed+ 10000;
    auto cseed = seed + 10001;
    auto eseed = seed+ 100000;

    auto weights = satgirgs::generateWeights(n);
    auto nc_positions = satgirgs::generatePositions(n, 2, ncseed);
    auto nc_nodes = satgirgs::convertToNodes(nc_positions, weights);
    auto c_positions = satgirgs::generatePositions(m, 2, cseed);
    std::vector<double> c_pseudoweights(m, 1); // clause nodes all have weight 1 in the model
    auto c_nodes = satgirgs::convertToNodes(c_positions, c_pseudoweights, nc_nodes.size());

   auto edges = satgirgs::generateEdges(c_nodes, nc_nodes, k, t, eseed, true);
   auto edgeCount = edges.size();
   auto variablesPerClause = edgeCount / m;

    cout << n << ','
         << m << ','
         << k << ','
         << t << ','
         << threads << ','
         << seed << ','
         << plot << ','
         << edgeCount << ','
         << variablesPerClause << '\n';
}


int main(int argc, char* argv[]) {

    // cout << "dimension,n,avgDeg,alpha,ple,threads,seed,plot,TimeWeights,TimePositions,TimeBinary,TimePre,TimeEdges,TimeTotal,GenNumEdge,GenAvgDeg\n";

    int seed = 0;

    auto n = 1<<15;
    auto m = 1<<10;
    auto t = 0.5;
    auto ple = 2.5;
    auto k=10;
    auto threads = 1;
    auto reps = 10;

    for(int rep=0; rep<reps; ++rep) {
        clog << "rep " << rep << endl;

        clog << "growing n" << endl;
        for(int i = 1<<10; i<= (1<<15); i <<= 1) {
            clog << i << endl;
            for(auto deg : {10,15, 20})
		if(deg*4 < i) 
		    measure(i, m, deg, t, threads, ++seed, 0);
        }


        clog << "growing deg" << endl;
        for(int i = 2; i<= 1<<6; i <<= 1) {
            clog << i << endl;
            measure(n, m, i, t, threads, ++seed, 1);
        }

        clog << "shrinking t" << endl;
        for(auto i : {5.0, 4.0, 3.0, 2.0, 1.0, 0.5, 0.0}) {
            clog << i << endl;
            measure(n, m, k, i, threads, ++seed, 2);
        }
    }

    return 0;
}
