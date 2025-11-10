#pragma once

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <map>

#include "lemon/list_graph.h"

namespace vgpart {
    using Graph = lemon::ListGraph;

    using EdgeKey = std::pair<uint32_t,uint32_t>;

    // 规范化边键（小id在前）
    inline EdgeKey norm_edge(uint32_t a, uint32_t b) {
      if (a < b) return {a, b};
      if (a > b) return {b, a};
      return {a, b};
    }

    struct SubsetPruneResult {
      std::vector<uint32_t> kept_nodes;
      std::vector<uint32_t> removed_nodes;
      std::vector<EdgeKey> bridges;
      bool has_bridge = false;
      int k_used = 2;
      int min_edge_weight = 0;

      std::vector<std::vector<uint32_t>> components;
      double lcc_ratio = 1.0;
    };

    struct SubsetPack {
      // 当前子集的“保留节点（LCC）”
      std::vector<uint32_t> core_nodes;
      // 当前子集的“孤分量”（不在 LCC 的各个分量）
      std::vector<std::vector<uint32_t>> orphan_components;
    };

    struct EdgeFusionOptions {
        double alpha = 0.6;
        double beta = 0.4;
        double sigma = -1.0;
        double average_ends = true;
        double min_edge_weight = 1e-6;
    };

    struct PruneCoresOptions {
        int k_core = 2;
        int min_edge_weight = 0;
        bool detect_bridges = false;
    };

    struct ReassignOptions {
        int k_core = 2;
        int min_edge_weight = 0;
        double keep_ratio_thresh = 1;
        long long score_floor = 0;
    };

    struct OverlapOptions {
      int K = 20;
      int k_core = 2;
      int min_edge_weight = 0;
      int top_blocks = 2;
    };
    }