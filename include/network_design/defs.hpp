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
    Graph graph;

    /// Função de custo das arestas.
    Graph::EdgeMap<double> edge_cost;

    /// Função de peso dos vértices.
    Graph::NodeMap<double> node_weight;

    /// Capacidade máxima das partições.
    double capacity;

    /// Fator multiplicador de custo das arestas no circuito (Gama).
    double circuit_cost_factor;

    instance_t() = delete;

    // Construtor de cópia desabilitado porque copiar a instância quebra os
    // mapas apontando para o grafo.
    instance_t(const instance_t&) = delete;

    instance_t(int nnodes);
};

/**
 * Estrutura para os dados de uma solução do problema.
 */
struct solution_t {
    const instance_t& instance;

    /// Partição dos nós do grafo.
    Graph::NodeMap<int> partition;

    /// Nós do circuito, em ordem.
    std::vector<Graph::Node> circuit_nodes;

    solution_t() = delete;
    solution_t(const instance_t&);
    solution_t(const solution_t& other);

    /// Gera os arcos do circuito na solução.
    std::vector<Graph::Arc> circuit_arcs() const;

    /// Gera os arcos das estrelas na solução.
    std::vector<Graph::Arc> star_arcs() const;

    /// Computa o custo da solução
    double cost() const;
};

}; // namespace network_design

#endif // _NETWORK_DESIGN_DEFS_HPP
