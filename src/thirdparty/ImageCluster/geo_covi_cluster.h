#pragma once

#include "graph_types.h"
#include "lemon_util.h"

#include "colmap/scene/scene_clustering.h"

namespace vgpart {
    struct GeoCoviClusterOptions {
      EdgeFusionOptions fusion_opts;
      PruneCoresOptions prune_opts;
      ReassignOptions reassign_opts;
      OverlapOptions overlap_opts;
    };
    
    // TODO:处理POSE == COVI场景 √
    // TODO:处理POSE < COVI场景
    // TODO:处理POSE = 0, COVI > 0场景
    class GeoCoviCluster {
     public:
      bool FuseGeoCoviEdgeWeight(
          const std::map<EdgeKey, int>& covi_edges_weight,
          const std::map<EdgeKey, int>& geo_edges_weight,
          const EdgeFusionOptions& fusion_opts,
          std::map<EdgeKey, int>& fusion_edges_weight);

      bool RunMetisSplit(const std::map<EdgeKey, int>& edges_weight,
                         const uint32_t cluster_size,
                         const uint32_t overlap_size,
                         std::vector<std::vector<uint32_t>>& output_clusters);

      bool RunPruneCores(
          const std::vector<std::vector<uint32_t>>& output_clusters,
          const std::map<EdgeKey, int>& edges_weight,
          const PruneCoresOptions prune_opts,
          std::vector<SubsetPack>& packs);

      bool RunReassignment(const std::map<EdgeKey, int>& edges_weight,
                           const ReassignOptions& reassign_opts,
                           std::vector<SubsetPack>& packs,
                           std::vector<std::vector<uint32_t>>& output_clusters);

      bool RunComputeOverlaps(
          const std::map<EdgeKey, int>& edges_weight,
          const OverlapOptions& overlap_opt,
          std::vector<std::vector<uint32_t>>& output_clusters);

      bool RunAll(const std::map<EdgeKey, int>& covi_edges_weight,
                  const std::map<EdgeKey, int>& geo_edges_weight,
                  const uint32_t cluster_size,
                  const uint32_t overlap_size,
                  GeoCoviClusterOptions& opts,
                  std::vector<std::vector<uint32_t>>& output_clusters);

      bool FillCluster(const std::vector<std::vector<uint32_t>>& output_clusters,
          colmap::SceneClustering::Cluster& root_cluster);

      bool SavePLY_XYZ(const std::string& path,
                       const std::vector<Eigen::Vector3d>& pts);

     private:
      // ============== 图分割 ===============
      bool PrepareMetisData(const std::map<EdgeKey, int>& edges_weight,
                            std::vector<EdgeKey>& vec_image_pairs,
                            std::vector<int>& pair_weights);

      bool FiltEdges(const std::map<EdgeKey, int>& input,
                     std::map<EdgeKey, int>& output);

      bool MetisSplit(const std::vector<EdgeKey>& image_pairs,
                      const std::vector<int> pair_weights,
                      uint32_t cluster_size,
                      uint32_t overlap_size,
                      std::vector<std::vector<uint32_t>>& output_clusters);
    };


}