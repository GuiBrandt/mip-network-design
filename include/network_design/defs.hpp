#ifndef _NETWORK_DESIGN_DEFS_HPP
#define _NETWORK_DESIGN_DEFS_HPP

#include <cstdint>

#include <lemon/full_graph.h>

namespace network_design {

using Graph = lemon::FullGraph;

/**
 * Estrutura para os dados de uma instância do problema.
 */
struct instance_t {
    /// Grafo da instância.
    const Graph& graph;

    /// Função de custo das arestas.
    Graph::EdgeMap<uint64_t> edge_cost;

    /// Função de peso dos vértices.
    Graph::NodeMap<uint64_t> node_weight;

    /// Capacidade máxima das partições.
    uint64_t capacity;

    /// Fator multiplicador de custo das arestas no circuito (Gama).
    uint64_t circuit_cost_factor;

    instance_t() = delete;
    instance_t(const Graph&);
    instance_t(const instance_t&);
};

/**
 * Estrutura para os dados de uma solução do problema.
 */
struct solution_t {
    const instance_t& instance;

    /// Partição dos nós do grafo.
    Graph::NodeMap<int> partition;

    /// Vetor de representantes das partições (na ordem numérica das partições)
    std::vector<Graph::Node> partition_repr;

    /// Nós do circuito, em ordem.
    std::vector<Graph::Node> circuit_nodes;

    solution_t() = delete;
    solution_t(const instance_t&);

    solution_t(const solution_t& other);
    solution_t(const solution_t&& other);

    /// Gera as arestas do circuito na solução.
    std::vector<Graph::Edge> circuit_edges() const;

    /// Gera as arestas das estrelas na solução.
    std::vector<Graph::Edge> star_edges() const;

    /// Computa o custo da solução
    uint64_t cost() const;
};

}; // namespace network_design

#endif // _NETWORK_DESIGN_DEFS_HPP
