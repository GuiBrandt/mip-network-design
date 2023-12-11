#include <iostream>

#include <mylib/myutils.h>
#include <network_design.hpp>

int main(int argc, char* argv[]) {
    using formulation_t = network_design::mip::branch_and_cut::formulation_t;

    if (argc != 2) {
        std::cout << std::endl
                  << "Usage: " << argv[0] << " <graph_filename>" << std::endl
                  << std::endl
                  << "Example:" << std::endl
                  << "\t" << argv[0] << " "
                  << getpath(argv[0]) +
                         "../instances/mo420_network_design_1_10_100_20_80_1"
                  << std::endl
                  << "\t" << argv[0] << " "
                  << getpath(argv[0]) +
                         "../instances/mo420_network_design_1_20_100_20_80_2"
                  << std::endl
                  << std::endl;
        std::exit(1);
    } else if (!FileExists(argv[1])) {
        std::cout << "File " << argv[1] << " does not exist." << endl;
        std::exit(1);
    }

    auto instance = network_design::read_instance(argv[1]);

    std::random_device rd;
    int seed = rd();
    std::cout << "Seed: " << seed << std::endl;

    std::minstd_rand rng(seed);
    std::uniform_int_distribution<int> grb_seed_dist(0, GRB_MAXINT);

    auto greedy_solution = network_design::greedy_heuristic(*instance);
    auto cutoff = greedy_solution.cost();
    std::cout << "Heuristic cost: " << cutoff << std::endl;

    GRBEnv env;
    env.set(GRB_IntParam_Seed, grb_seed_dist(rng));
    env.set(GRB_DoubleParam_Cutoff, cutoff);
    formulation_t formulation(*instance, env);

    auto solution = formulation.solve();
    network_design::view(solution);

    return 0;
}
