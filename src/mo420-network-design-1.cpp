#include <chrono>
#include <iostream>

#include <mylib/myutils.h>
#include <network_design.hpp>

int main(int argc, char* argv[]) {
    using formulation_t = network_design::mip::polynomial::formulation_t;

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
        std::cout << "File " << argv[1] << " does not exist." << std::endl;
        std::exit(1);
    }

    auto instance_file = std::string(argv[1]);
    auto instance_name =
        instance_file.substr(instance_file.find_last_of("/\\") + 1);
    auto instance = network_design::read_instance(instance_file);

    std::random_device rd;
    int seed = rd();
    std::cout << "Seed: " << seed << std::endl;

    std::minstd_rand rng(seed);
    std::uniform_int_distribution<int> grb_seed_dist(0, GRB_MAXINT);

    auto start_time = std::chrono::steady_clock::now();

    auto greedy_solution = network_design::greedy_heuristic(*instance);
    auto cutoff = greedy_solution.cost();
    std::cout << "Heuristic cost: " << cutoff << std::endl;
    std::cout << "Heuristic time: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(
                     std::chrono::steady_clock::now() - start_time)
                     .count()
              << "ms" << std::endl;

    GRBEnv env;
    env.set(GRB_IntParam_Seed, grb_seed_dist(rng));
    env.set(GRB_DoubleParam_Cutoff, cutoff);

    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now() - start_time);
    env.set(GRB_DoubleParam_TimeLimit, 3600 - elapsed.count() / 1000.0);

    formulation_t formulation(*instance, env);

    auto solution = formulation.solve();
    network_design::view(solution, instance_name);

    std::cout << "Total time: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(
                     std::chrono::steady_clock::now() - start_time)
                     .count()
              << "ms" << std::endl;

    return 0;
}
