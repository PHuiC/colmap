#include "hierarchical_cluster.h"

bool vgpart::HierarchicalCluster::RunPartition(
    const std::map<EdgeKey, int>& covi_edges_weight,
    const std::map<EdgeKey, int>& geo_edges_weight,
    HierarchicalOptions& opts,
    const int level,
    colmap::SceneClustering::Cluster* root_cluster) {

    root_cluster->level = level;

    if (covi_edges_weight.empty() ||
      root_cluster->image_ids.size() <=
                           static_cast<size_t>(opts.leaf_max_num_images)) {
        return true;
    }
    std::vector<std::vector<uint32_t>> geo_covi_result;
    geo_covi_cluster_.RunAll(covi_edges_weight,
                             geo_edges_weight,
                             opts.branching,
                             opts.image_overlap,
                             opts.geo_covi_options,
                             geo_covi_result);

    geo_covi_cluster_.FillCluster(geo_covi_result, *root_cluster);

    std::vector<std::map<EdgeKey, int>> child_covi_edges_weight(opts.branching);
    std::vector<std::map<EdgeKey, int>> child_geo_edges_weight(opts.branching);
    
    for (int i = 0; i < opts.branching; ++i) {
      std::set<uint32_t> set_cluster(geo_covi_result[i].begin(),
                                     geo_covi_result[i].end());

      for (auto& kv : covi_edges_weight) {
        auto& pair = kv.first;
        auto& weight = kv.second;
        if (set_cluster.count(pair.first) > 0 && 
            set_cluster.count(pair.second) > 0) {
          child_covi_edges_weight[i].insert({pair, weight});
        }
      }

      for (auto& kv : geo_edges_weight) {
        auto& pair = kv.first;
        auto& weight = kv.second;
        if (set_cluster.count(pair.first) > 0 &&
            set_cluster.count(pair.second) > 0) {
          child_geo_edges_weight[i].insert({pair, weight});
        }
      }
    }

    for (int i = 0; i < opts.branching; ++i) {
      if (root_cluster->child_clusters[i].image_ids.empty() ||
          root_cluster->child_clusters[i].image_ids.size() ==
              root_cluster->image_ids.size()) {
        continue;
      }

      RunPartition(child_covi_edges_weight[i],
                   child_geo_edges_weight[i],
                   opts,
                   level + 1,
                   &root_cluster->child_clusters[i]);
    }

    // Remove empty clusters.
    root_cluster->child_clusters.erase(
        std::remove_if(root_cluster->child_clusters.begin(),
                       root_cluster->child_clusters.end(),
            [](const colmap::SceneClustering::Cluster& childCluster) {
                         return childCluster.image_ids.empty();
                       }),
        root_cluster->child_clusters.end());

    // If the child cluster is the same as the current cluster, it is redundant
    // and we can remove it.
    if (root_cluster->child_clusters.size() == 1 &&
        root_cluster->image_ids.size() ==
            root_cluster->child_clusters[0].image_ids.size()) {
      root_cluster->child_clusters = {};
    }


    return true;
}
