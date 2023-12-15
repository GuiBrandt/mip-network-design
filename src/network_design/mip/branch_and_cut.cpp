#include "network_design/mip/branch_and_cut.hpp"

#include <limits>

#include <lemon/gomory_hu.h>
#include <lemon/list_graph.h>
#include <lemon/unionfind.h>

namespace network_design {
namespace mip {
namespace branch_and_cut {

formulation_t::formulation_t(const instance_t& instance, const GRBEnv& env)
    : formulation_base(instance, env) {
    model.set(GRB_IntParam_LazyConstraints, 1);
    model.setCallback(this);
}

static char constr_name[64];

void formulation_t::find_violated_integer_cuts() {
    const auto& G = instance.graph;

    Graph::EdgeMap<double> capacity(G);
    for (Graph::EdgeIt e(G); e != lemon::INVALID; ++e) {
        Graph::Node u = G.u(e), v = G.v(e);
        Graph::Arc a = G.arc(u, v), b = G.arc(v, u);
        capacity[e] = getSolution(vars.circuit_arc[a]) +
                      2 * getSolution(vars.star_arc[a]) +
                      getSolution(vars.circuit_arc[b]) +
                      2 * getSolution(vars.star_arc[b]);
    }

    Graph::NodeMap<int> aux_map(G);
    lemon::UnionFind<Graph::NodeMap<int>> components(aux_map);
    for (Graph::NodeIt v(G); v != lemon::INVALID; ++v) {
        components.insert(v);
    }
    for (Graph::EdgeIt e(G); e != lemon::INVALID; ++e) {
        if (capacity[e] > 1 - 1e-4) {
            components.join(G.u(e), G.v(e));
        }
    }

    int num_components = 0;
    std::vector<int> component_index(G.nodeNum(), -1);
    for (Graph::NodeIt v(G); v != lemon::INVALID; ++v) {
        int component = components.find(v);
        if (component_index[component] < 0) {
            component_index[component] = num_components++;
        }
    }

    if (num_components == 1) {
        return;
    }

    for (int i = 0; i < num_components; i++) {
        GRBLinExpr out_expr, in_expr;
        for (Graph::ArcIt a(G); a != lemon::INVALID; ++a) {
            if (component_index[components.find(G.source(a))] == i &&
                component_index[components.find(G.target(a))] != i) {
                out_expr += vars.circuit_arc[a] + vars.star_arc[a] +
                            vars.star_arc[G.oppositeArc(a)];
            }
        }
        addLazy(out_expr >= 1);
    }
}

void formulation_t::find_violated_fractional_cuts() {
    const auto& G = instance.graph;

    Graph::EdgeMap<double> capacity(G);
    for (Graph::EdgeIt e(G); e != lemon::INVALID; ++e) {
        Graph::Node u = G.u(e), v = G.v(e);
        Graph::Arc a = G.arc(u, v), b = G.arc(v, u);
        capacity[e] =
            getNodeRel(vars.circuit_arc[a]) + 2 * getNodeRel(vars.star_arc[a]) +
            getNodeRel(vars.circuit_arc[b]) + 2 * getNodeRel(vars.star_arc[b]);
    }

    std::vector<Graph::Edge> frac_edges, one_edges;
    for (Graph::EdgeIt e(G); e != lemon::INVALID; ++e) {
        if (capacity[e] > 1 - 1e-4) {
            one_edges.push_back(e);
        } else if (capacity[e] > 1e-4) {
            frac_edges.push_back(e);
        }
    }

    Graph::NodeMap<int> aux_map(G);
    lemon::UnionFind<Graph::NodeMap<int>> components(aux_map);
    for (Graph::NodeIt v(G); v != lemon::INVALID; ++v) {
        components.insert(v);
    }
    for (auto e : one_edges) {
        components.join(G.u(e), G.v(e));
    }

    int num_components = 0;
    std::vector<int> component_index(G.nodeNum(), -1);
    for (Graph::NodeIt v(G); v != lemon::INVALID; ++v) {
        int component = components.find(v);
        if (component_index[component] < 0) {
            component_index[component] = num_components++;
        }
    }

    using LGraph = lemon::ListGraph;

    LGraph H;

    std::vector<LGraph::Node> H_index(num_components);
    for (int i = 0; i < num_components; i++) {
        H_index[i] = H.addNode();
    }

    LGraph::EdgeMap<double> H_capacity(H);
    for (auto e : frac_edges) {
        auto u = G.u(e), v = G.v(e);
        auto hu = H_index[component_index[components.find(u)]],
             hv = H_index[component_index[components.find(v)]];
        auto he = H.addEdge(hu, hv);
        H_capacity[he] = capacity[e];
    }

    lemon::GomoryHu<LGraph, LGraph::EdgeMap<double>> gomory_hu(H, H_capacity);
    gomory_hu.run();

    LGraph::NodeMap<bool> cutmap(H, false);
    for (LGraph::NodeIt u(H); u != lemon::INVALID; ++u) {
        if (gomory_hu.predNode(u) == lemon::INVALID ||
            gomory_hu.predValue(u) > 2.0 - 1e-5) {
            continue;
        }

        gomory_hu.minCutMap(u, gomory_hu.predNode(u), cutmap);

        GRBLinExpr out_expr;
        for (Graph::ArcIt a(G); a != lemon::INVALID; ++a) {
            auto s = G.source(a), t = G.target(a);
            auto hs = H_index[component_index[components.find(s)]],
                 ht = H_index[component_index[components.find(t)]];
            if (cutmap[hs] && !cutmap[ht]) {
                out_expr += vars.circuit_arc[a] + vars.star_arc[a] +
                            vars.star_arc[G.oppositeArc(a)];
            }
        }

        addCut(out_expr >= 1);
    }
}

void formulation_t::find_violated_blossom() {
    const auto& G = instance.graph;

    Graph::EdgeMap<double> capacity(G);
    for (Graph::EdgeIt e(G); e != lemon::INVALID; ++e) {
        Graph::Node u = G.u(e), v = G.v(e);
        Graph::Arc a = G.arc(u, v), b = G.arc(v, u);
        capacity[e] =
            getNodeRel(vars.circuit_arc[a]) + 2 * getNodeRel(vars.star_arc[a]) +
            getNodeRel(vars.circuit_arc[b]) + 2 * getNodeRel(vars.star_arc[b]);
    }

    lemon::GomoryHu<Graph, Graph::EdgeMap<double>> gomory_hu(G, capacity);
    gomory_hu.run();

    Graph::NodeMap<int> in_degree(G);
    Graph::NodeMap<int> cut_size(G);

    std::list<Graph::Node> sources, topological_sort;
    for (Graph::NodeIt v(G); v != lemon::INVALID; ++v) {
        in_degree[v] = 0;
        cut_size[v] = 1;
    }

    for (Graph::NodeIt u(G); u != lemon::INVALID; ++u) {
        auto v = gomory_hu.predNode(u);
        if (v == lemon::INVALID) {
            continue;
        }
        in_degree[v]++;
    }

    for (Graph::NodeIt v(G); v != lemon::INVALID; ++v) {
        if (in_degree[v] == 0) {
            sources.push_back(v);
        }
    }

    while (!sources.empty()) {
        auto u = sources.front();
        sources.pop_front();
        topological_sort.push_back(u);
        auto v = gomory_hu.predNode(u);
        if (v != lemon::INVALID) {
            in_degree[v]--;
            if (in_degree[v] == 0) {
                sources.push_back(v);
            }
        }
    }

    while (!topological_sort.empty()) {
        auto u = topological_sort.front();
        topological_sort.pop_front();
        auto v = gomory_hu.predNode(u);
        if (v != lemon::INVALID) {
            cut_size[v] = cut_size[v] + cut_size[u];
        }
    }

    for (Graph::NodeIt u(G); u != lemon::INVALID; ++u) {
        GRBLinExpr expr;
        if (cut_size[u] % 2 == 0 || gomory_hu.predValue(u) > 1.0 - 1e-5) {
            continue;
        }

        for (decltype(gomory_hu)::MinCutEdgeIt e(gomory_hu, u,
                                                 gomory_hu.predNode(u));
             e != lemon::INVALID; ++e) {
            Graph::Arc a = G.arc(G.u(e), G.v(e)), b = G.oppositeArc(a);
            expr += vars.circuit_arc[a] + vars.circuit_arc[b] +
                    2 * (vars.star_arc[a] + vars.star_arc[b]);
        }

        addCut(expr >= 1.0);
    }
}

void formulation_t::callback() {
    try {
        switch (where) {
        case GRB_CB_MIPSOL:
            find_violated_integer_cuts();
            break;

        case GRB_CB_MIPNODE:
            if (getIntInfo(GRB_CB_MIPNODE_STATUS) == GRB_OPTIMAL) {
                // find_violated_blossom();
                find_violated_fractional_cuts();
            }
            break;

        default:
            break;
        }
    } catch (const GRBException& ex) {
        std::cerr << "ERROR: " << ex.getMessage() << " (" << ex.getErrorCode()
                  << ")" << std::endl;
    }
}

} // namespace branch_and_cut
} // namespace mip
} // namespace network_design
