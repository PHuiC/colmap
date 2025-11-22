#pragma once

#include "graph_types.h"
#include "lemon_util.h"

#include "geo_covi_cluster.h"

#include "colmap/scene/scene_clustering.h"


namespace vgpart {
struct HierarchicalOptions {
  GeoCoviClusterOptions geo_covi_options;

  int branching = 2;

  int image_overlap = 50;

  int leaf_max_num_images = 1200;
};

class HierarchicalCluster {
 public:
  bool RunPartition(const std::map<EdgeKey, int>& covi_edges_weight,
           const std::map<EdgeKey, int>& geo_edges_weight,
           HierarchicalOptions& opts,
           const int level,
           colmap::SceneClustering::Cluster* root_cluster);

 private:
  GeoCoviCluster geo_covi_cluster_;
};
}  // namespace vgpart