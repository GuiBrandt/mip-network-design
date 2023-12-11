#include <network_design.hpp>

int main(int argc, char* argv[]) {
    using formulation_t = network_design::mip::polynomial::formulation_t;

    std::random_device rd;
    int seed = rd();
    std::cout << "seed: " << seed << std::endl;

    std::minstd_rand rng(seed);
    std::uniform_int_distribution<int> grb_seed_dist(0, GRB_MAXINT);

    network_design::Graph G(50);
    auto instance = network_design::random_instance(G, seed);

    auto greedy_solution = network_design::greedy_heuristic(instance);
    auto cutoff = greedy_solution.cost();
    std::cout << "Heuristic cost: " << cutoff << std::endl;

    GRBEnv env;
    env.set(GRB_IntParam_Seed, grb_seed_dist(rng));
    env.set(GRB_DoubleParam_Cutoff, cutoff);
    formulation_t formulation(instance, env);

    auto solution = formulation.solve();
    network_design::view(solution);

    return 0;
}
