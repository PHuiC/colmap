#pragma once

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <map>
#include <algorithm>
#include <cstdint>

#include "graph_types.h"

namespace vgpart {

    SubsetPruneResult PruneSubsetByKcoreAndBridges(
        const std::vector<uint32_t>& subset_nodes,
        const std::map<EdgeKey, int>& edges_weight,
        int k_core = 2,
        int min_edge_weight = 0,
        bool detect_bridges = true);

    void BuildGraphForSubset(
        const std::vector<uint32_t>& subset_nodes,
        const std::map<EdgeKey, int>& edges_weight,
        int min_edge_weight,
        Graph& g,
        std::unordered_map<uint32_t, Graph::Node>& id2node,
        lemon::ListGraph::NodeMap<uint32_t>& node_id_map,
        std::vector<Graph::Edge>& edge_list,
        std::vector<EdgeKey>& edge_keys);

    void KCorePrune(Graph& g,
                      const lemon::ListGraph::NodeMap<uint32_t>& node_id_map,
                      int k,
                      std::vector<uint32_t>& kept,
                      std::vector<uint32_t>& removed);

    void FindBridges(const Graph& g,
                      const lemon::ListGraph::NodeMap<uint32_t>& node_id_map,
                      const std::vector<Graph::Edge>& edge_list,
                      const std::vector<EdgeKey>& edge_keys,
                      std::vector<EdgeKey>& bridges_out);

    std::vector<std::vector<uint32_t>> ExtractConnectedComponents(
        const std::vector<uint32_t>& kept_nodes,
        const std::map<EdgeKey, int>& edges_weight,
        int min_edge_weight);

    long long CrossScoreComponentToSubset(
        const std::vector<uint32_t>& comp,
        const std::vector<uint32_t>&
            target_subset_nodes,  // 用剪枝后 LCC 的节点
        const std::map<EdgeKey, int>& edges_weight,
        int min_edge_weight);

    bool KcoreFeasibleAfterMerge(
        const std::vector<uint32_t>& comp,
        const std::vector<uint32_t>& target_subset_nodes,
        const std::map<EdgeKey, int>& edges_weight,
        int k_core,
        int min_edge_weight,
        double keep_ratio_thresh);

    bool ReassignComponentsGreedy(std::vector<SubsetPack>& subsets,  // in/out
                                    const std::map<EdgeKey, int>& edges_weight,
                                    const ReassignOptions& opt);

    bool ReassignPoolImagewise(
        std::vector<SubsetPack>& subsets,  // in/out
        const std::map<EdgeKey, int>& edges_weight,
        const ReassignOptions& opt);

    int EffectiveDegreeToSet(uint32_t u,
                                const std::unordered_set<uint32_t>& T,
                                const std::map<EdgeKey, int>& edges_weight,
                                int min_edge_weight);

    long long CrossScoreToSet(uint32_t u,
                                 const std::unordered_set<uint32_t>& T,
                                 const std::map<EdgeKey, int>& edges_weight,
                                 int min_edge_weight);

    std::vector<uint32_t> SelectOverlapFromSrcToDst(
        const std::vector<uint32_t>& A_core,
        const std::vector<uint32_t>& B_core,
        const std::map<EdgeKey, int>& edges_weight,
        const OverlapOptions& p);

    std::vector<std::vector<std::vector<uint32_t>>> ComputePairwiseOverlaps(
        const std::vector<std::vector<uint32_t>>& cores,
        const std::map<EdgeKey, int>& edges_weight,
        const OverlapOptions& p);
    }