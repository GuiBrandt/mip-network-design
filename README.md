# PLI - Projetos (Network Design)

O projeto está organizado da seguinte forma:
- `CMakeLists.txt` - Arquivo com instruções de compilação para o CMake.
- `Install-Lemon-in-subfolder.bash` - Script para instalação da biblioteca LEMON.
- `include/` - Pasta com os arquivos de cabeçalho.
  - `mylib/` - Cabeçalhos da biblioteca auxiliar.
  - `network_design/` - Cabeçalhos dos algoritmos relacionados ao problema.
    - `mip/` - Cabeçalhos dos algoritmos de programação inteira.
      - `base.hpp` - Núcleo comum das formulações nas duas partes do projeto.
      - `branch_and_cut.hpp` - Formulação por cortes.
      - `polynomial.hpp` - Formulação compacta (de tamanho polinomial).
    - `defs.hpp` - Definição das estruturas de instância e solução do problema.
    - `heuristics.hpp` - Definição da heurística gulosa.
    - `utils.hpp` - Outras funções (e.g. leitura de instância, visualizar PDF, etc.)
- `src/` - Pasta com os arquivos de código C++.
  - `mylib/` - Código fonte da biblioteca auxiliar.
  - `network_design/` - Código fonte dos algoritmos relacionados ao problema. A estrutura corresponde à da pasta `include/network_design`.
  - `mo420-network-design-1.cpp` - Arquivo principal do programa da parte 1.
  - `mo420-network-design-2.cpp` - Arquivo principal do programa da parte 2.
- `instances/` - Instâncias do problema utilizadas para teste.
- `report/` - Relatórios.
  - `polynomial/main.pdf` - Relatório da parte 1.
  - `branch-and-cut/main.pdf` - Relatório da parte 2.

O projeto pode ser compilado com uma versão compatível do CMake (3.22+) e um compilador com suporte a C++ 11. Testado com CMake 3.28 e Clang 16.0.6.

## Parte 1 - Formulação compacta

$$
\begin{align}
    (\mathcal{P})\quad
    &\min \sum_{(u, v) \in A} \gamma c_{uv} x_{uv} + \sum_{(u, v) \in A} c_{uv} y_{uv}\nonumber\\
    &\sum_{v \in V} r_{v} \geq 3\\
    &\sum_{v \in \delta^-(u)} y_{vu} = 1 - r_u              &\forall u \in V\\
    &\sum_{v \in \delta^+(u)} d_v y_{uv} \leq r_u (B - d_u) &\forall u \in V\\
    &\sum_{v \in \delta^+(u)} x_{uv} = r_u                  &\forall u \in V\\
    &\sum_{v \in \delta^-(u)} x_{vu} = r_u                  &\forall u \in V\\
    &t_u - t_v + |V|(x_{uv} - y_{vv_0}) \leq |V| - 1&\forall (u, v) \in A, v \neq v_0\\
    &x_{uv} \in \{0, 1\}&\forall (u, v) \in A\nonumber\\
    &y_{uv} \in \{0, 1\}&\forall (u, v) \in A\nonumber\\
    &r_{v} \in \{0, 1\}&\forall v \in V\nonumber\\
    &t_{v} \in \mathbb{R}_{\geq 0}&\forall v \in V\nonumber
\end{align}
$$

## Parte 2 - Branch-and-Cut

$$
\begin{align}
    (\mathcal{P})\quad
    &\min \sum_{(u, v) \in A} \gamma c_{uv} x_{uv} + \sum_{(u, v) \in A} c_{uv} y_{uv}\nonumber\\
    &\sum_{v \in V} r_{v} \geq 3
        \\
    &\sum_{v \in \delta^-(u)} y_{vu} = 1 - r_u
        &\forall u \in V\\
    &\sum_{v \in \delta^+(u)} d_v y_{uv} \leq r_u (B - d_u)
        &\forall u \in V\\
    &\sum_{v \in \delta^+(u)} x_{uv} = r_u
        &\forall u \in V\\
    &\sum_{v \in \delta^-(u)} x_{vu} = r_u
        &\forall u \in V\\
    &\sum_{(u, v) \in \delta^+(S)} x_{uv} + y_{uv} + y_{vu} \geq 1
        &\forall S \subsetneq V\\
    &x_{uv} \in \{0, 1\}&\forall (u, v) \in A\nonumber\\
    &y_{uv} \in \{0, 1\}&\forall (u, v) \in A\nonumber\\
    &r_{v} \in \{0, 1\}&\forall v \in V\nonumber
\end{align}
$$
