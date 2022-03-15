
#include <iostream>
#include <chrono>
#include <map>
#include <string>
#include <cstring>
#include <fstream>
#include <algorithm>

#include <omp.h>

#include <girgs/girgs-version.h>
#include <satgirgs/Generator.h>
#include <satgirgs/BitManipulation.h>


using namespace std;
using namespace chrono;


map<string, string> parseArgs(int argc, char** argv) {
    map<string, string> params;
    for (int i = 1; i < argc; i++) {
        // Get current and next argument
        if (argv[i][0] != '-')
            continue;
        std::string arg = argv[i] + 1; // +1 to skip the -
        // advance one additional position if next is used
        std::string next = (i + 1 < argc ? argv[i++ + 1] : ""); 
        params[std::move(arg)] = std::move(next);
    }
    return params;
}


template<typename T>
void logParam(T value, const string& name) {
    cout << "\t" << name << "\t=\t" << value << '\n';
}

template<typename T>
void rangeCheck(T value, T min, T max, string name, bool lex = false, bool hex = false) {
    if (value < min || value > max || (value == min && lex) || (value == max && hex)) {
        cerr << "ERROR: parameter " << name << " = " << value << " is not in range "
            << (lex ? "(" : "[") << min << "," << max << (hex ? ")" : "]") << '\n';
        exit(1);
    }
    logParam(value, name);
}



int main(int argc, char* argv[]) {

    // write help
    if (argc < 2 || 0 == strcmp(argv[1], "--help") || 0 == strcmp(argv[1], "-help")) {
        clog << "usage: ./gensatgirg\n"
            << "\t\t[-n anInt]          // number of vertices (non-clause points)   default 10000\n"
            << "\t\t[-m anInt]          // number of edges (clause points)          default 10000\n"
            << "\t\t[-ple aFloat]       // power law exponent       range (2,3]     default 2.5\n"
            << "\t\t[-alpha aFloat]     // model parameter          range (1,inf]   default infinity\n"
            << "\t\t[-wseed anInt]      // weight seed                              default 12\n"
            << "\t\t[-ncseed anInt]     // non-clause position seed                 default 130\n"
            << "\t\t[-cseed anInt]      // clause position seed                     default 130\n"
            << "\t\t[-threads anInt]    // number of threads to use                 default 1\n"
            << "\t\t[-file aString]     // file name for output (w/o ext)           default \"graph\"\n"
            << "\t\t[-dot 0|1]          // write result as dot (.dot)               default 0\n"
            << "\t\t[-edge 0|1]         // write result as edgelist (.txt)          default 0\n"
            << "\t\t[-debug 0|1]         // output debug graph                      default 0\n";
        return 0;
    }

    // write version
    if(argc > 1 && 0 == strcmp(argv[1], "--version")) {
        cout << "SAT-GIRGs command line interface.\n\n"
             << GIRGS_NAME_VERSION << '\n'
             << GIRGS_PROJECT_DESCRIPTION << '\n'
             << GIRGS_AUTHOR_ORGANIZATION << '\n'
             << GIRGS_AUTHOR_DOMAIN << '\n'
             << GIRGS_AUTHOR_MAINTAINER << '\n';
        return 0;
    }

    // read params
    auto params = parseArgs(argc, argv);
    auto n      = !params["n"    ].empty()  ? stoi(params["n"    ]) : 10000;
    auto m      = !params["m"    ].empty()  ? stoi(params["m"    ]) : 10000; // TODO find sensible default and change in usage above
    auto ple    = !params["ple"  ].empty()  ? stod(params["ple"  ]) : 2.5;
    auto alpha  = !params["alpha"].empty()  ? stod(params["alpha"]) : std::numeric_limits<double>::infinity();
    auto wseed  = !params["wseed"].empty()  ? stoi(params["wseed"]) : 12;
    auto ncseed = !params["ncseed"].empty() ? stoi(params["ncseed"]): 130;
    auto cseed  = !params["cseed"].empty()  ? stoi(params["cseed"]) : 130;
    auto threads= !params["threads"].empty()? stoi(params["threads"]) : 1;
    auto file   = !params["file" ].empty()  ? params["file"] : "graph";
    auto dot    = params["dot" ] == "1";
    auto edge   = params["edge"] == "1";
    auto debug  = params["debug"] == "1";

    // log params and range checks
    cout << "using:\n";
    logParam(n, "n");
    logParam(m, "m");
    rangeCheck(ple, 2.0, 3.0, "ple", true, false);
    rangeCheck(alpha, 1.0, std::numeric_limits<double>::infinity(), "alpha", true);
    logParam(wseed, "wseed");
    logParam(ncseed, "ncseed");
    logParam(cseed, "cseed");
    rangeCheck(threads, 1, omp_get_max_threads(), "threads");
    omp_set_num_threads(threads);
    logParam(file, "file");
    logParam(dot, "dot");
    logParam(edge, "edge");
    logParam(debug, "debugMode");
    logParam(satgirgs::BitManipulation<1>::name(), "morton");
    cout << "\n";

    auto t1 = high_resolution_clock::now();


    cout << "generating weights ...\t\t" << flush;
    auto weights = satgirgs::generateWeights(n, ple, wseed, threads > 1);
    auto t2 = high_resolution_clock::now();
    cout << "done in " << duration_cast<milliseconds>(t2 - t1).count() << "ms\tlargest = ";
    cout << *max_element(weights.begin(), weights.end()) << endl;


    cout << "generating non-clause positions ...\t" << flush;
    auto nc_positions = satgirgs::generatePositions(n, 2, ncseed);
    auto nc_nodes = satgirgs::convertToNodes(nc_positions, weights);
    auto t3 = high_resolution_clock::now();
    cout << "done in " << duration_cast<milliseconds>(t3 - t2).count() << "ms" << endl;


    cout << "generating clause positions ...\t\t" << flush;
    auto c_positions = satgirgs::generatePositions(m, 2, cseed);
    std::vector<double> c_pseudoweights(m, 1);
    auto c_nodes = satgirgs::convertToNodes(c_positions, c_pseudoweights);
    auto t4 = high_resolution_clock::now();
    cout << "done in " << duration_cast<milliseconds>(t4 - t3).count() << "ms" << endl;

    cout << "sampling edges ...\t\t" << flush;
    auto edges = satgirgs::generateEdges(c_nodes, nc_nodes);
    if(debug) {
        auto debug_edges = satgirgs::generateEdges(c_nodes, nc_nodes, true);
    }
    auto t5 = high_resolution_clock::now();
    cout << "done in " << duration_cast<milliseconds>(t5 - t4).count() << "ms\tavg deg = " << edges.size()*2.0/n << endl;

    // TODO adapt saveDot code for new model
    if (dot) {
        cout << "writing .dot file ...\t\t" << flush;
        auto t6 = high_resolution_clock::now();
        satgirgs::saveDot(weights, nc_positions, edges, file+".dot");
        auto t7 = high_resolution_clock::now();
        cout << "done in " << duration_cast<milliseconds>(t7 - t6).count() << "ms" << endl;
    }

    if (edge) {
        cout << "writing edge list (.txt) ...\t" << flush;
        auto t6 = high_resolution_clock::now();
        ofstream f{file+".txt"};
        if(!f.is_open()) throw std::runtime_error{"Error: failed to open file \"" + file + ".txt\""};
        f << n << ' ' << edges.size() << "\n\n";
        for(auto& each : edges)
            f << each.first << ' ' << each.second << '\n';
        auto t7 = high_resolution_clock::now();
        cout << "done in " << duration_cast<milliseconds>(t7 - t6).count() << "ms" << endl;
    }

    return 0;
}
