#include "geo_covi_cluster.h"

#include <cassert>
#include <set>
#include <fstream>
#include <iomanip>

namespace vgpart {
    bool GeoCoviCluster::FuseGeoCoviEdgeWeight(
        const std::map<EdgeKey, int>& covi_edges_weight,
        const std::map<EdgeKey, int>& geo_edges_weight,
        const EdgeFusionOptions& fusion_opts,
        std::map<EdgeKey, int>& fusion_edges_weight) {
  
        double alpha = fusion_opts.alpha;
        double beta = fusion_opts.beta;
        double sigma = fusion_opts.sigma;
        double average_ends = fusion_opts.average_ends;
        double min_edge_weight = fusion_opts.min_edge_weight;

        std::map<EdgeKey, double> fused;

        // 1) 断言键集合一模一样（你已保证；Release 可去掉）
        assert(covi_edges_weight.size() == geo_edges_weight.size());
        {
          auto ita = covi_edges_weight.begin();
          auto itb = geo_edges_weight.begin();
          for (; ita != covi_edges_weight.end(); ++ita, ++itb) {
            assert(itb != geo_edges_weight.end() && itb->first == ita->first);
          }
        }

        // 2) sigma 自动估计（距离中位数）
        if (sigma <= 0.0) {
          std::vector<double> dists;
          dists.reserve(geo_edges_weight.size());
          for (auto& kv : geo_edges_weight)
            dists.push_back(static_cast<double>(kv.second));
          if (!dists.empty()) {
            size_t mid = dists.size() / 2;
            std::nth_element(dists.begin(), dists.begin() + mid, dists.end());
            double med = dists[mid];
            sigma = (med > 0.0) ? med : 1.0;
          } else {
            sigma = 1.0;
          }
        }

        // 3) 计算每个顶点的“本地最大 inliers”，用于归一化
        uint32_t max_vid = 0;
        for (auto& kv : covi_edges_weight)
          max_vid = std::max(max_vid, std::max(kv.first.first, kv.first.second));
        std::vector<int> local_max(max_vid + 1, 1);
        for (auto& kv : covi_edges_weight) {
          const EdgeKey& e = kv.first;
          int inl = kv.second;
          local_max[e.first] = std::max(local_max[e.first], inl);
          local_max[e.second] = std::max(local_max[e.second], inl);
        }

        // 4) 融合
        auto ita = covi_edges_weight.begin();
        auto itb = geo_edges_weight.begin();
        for (; ita != covi_edges_weight.end(); ++ita, ++itb) {
          const EdgeKey& e = ita->first;
          const int inl = ita->second;
          const double dij = static_cast<double>(itb->second);

          // s_ij：log1p 归一化，降低巨匹配边的支配力
          auto norm = [](int v, int mx) -> double {
            double a = std::log1p(static_cast<double>(v));
            double b = std::log1p(static_cast<double>(mx));
            return (b > 0.0) ? (a / b) : 0.0;
          };
          int max_i = (e.first < local_max.size()) ? local_max[e.first] : 1;
          int max_j = (e.second < local_max.size()) ? local_max[e.second] : 1;
          double si = norm(inl, std::max(1, max_i));
          double sj = norm(inl, std::max(1, max_j));
          double s_ij = average_ends ? 0.5 * (si + sj) : std::max(si, sj);

          // p_ij：高斯接近度
          double x = dij / sigma;
          double p_ij =
              (std::isfinite(x) && x < 1e6) ? std::exp(-0.5 * x * x) : 0.0;

          double w = alpha * s_ij + beta * p_ij;
          if (w >= min_edge_weight) fused.emplace(e, w);
        }

        double gamma = 1.0;
        for (auto& edge : fused) {
          const double S = 1000.0;  // 或 10000，看你想要的离散度
          int adjwgt = std::max(1, (int)std::lround(std::pow(edge.second, gamma) * S));
          fusion_edges_weight.insert({edge.first, adjwgt});
        }

    }

    bool GeoCoviCluster::RunMetisSplit(
        const std::map<EdgeKey, int>& edges_weight,
        const uint32_t cluster_size,
        const uint32_t overlap_size,
        std::vector<std::vector<uint32_t>>& output_clusters) {

      std::vector<EdgeKey> vec_image_pairs;
      std::vector<int> pair_weights;

      PrepareMetisData(edges_weight, vec_image_pairs, pair_weights);

      MetisSplit(vec_image_pairs,
                 pair_weights,
                 cluster_size,
                 overlap_size,
                 output_clusters);

      return true;
    }

    bool GeoCoviCluster::RunPruneCores(
        const std::vector<std::vector<uint32_t>>& output_clusters,
        const std::map<EdgeKey, int>& edges_weight,
        const PruneCoresOptions prune_opts,
        std::vector<SubsetPack>& packs) {

      packs.resize(output_clusters.size());
      for (size_t i = 0; i < output_clusters.size(); ++i) {
        auto R = PruneSubsetByKcoreAndBridges(output_clusters[i],
                                              edges_weight,
                                              prune_opts.k_core,
                                              prune_opts.min_edge_weight,
                                              prune_opts.detect_bridges);

        LOG(INFO) << "Cluster" << i
                  << " components size: " << R.components.size();
        LOG(INFO) << "Cluster" << i
                  << " removed nodes size: " << R.removed_nodes.size();

        if (!R.components.empty()) {
          packs[i].core_nodes = R.components[0];  // LCC
          for (size_t c = 1; c < R.components.size(); ++c)
            packs[i].orphan_components.push_back(std::move(R.components[c]));
        }

        if (!R.removed_nodes.empty()) {
          packs[i].orphan_components.push_back(std::move(R.removed_nodes));
        }
      }

      return true;
    }

    bool GeoCoviCluster::RunReassignment(
        const std::map<EdgeKey, int>& edges_weight,
        const ReassignOptions& reassign_opt,
        std::vector<SubsetPack>& packs,
        std::vector<std::vector<uint32_t>>& output_clusters) {
      bool ch1 = false;
      bool ch2 = true;
      while (ch1 || ch2) {
        // ch1 = ReassignComponentsGreedy(packs, min_edges_weight_filted, opt);
        LOG(INFO) << "Reassign...";
        ch2 = ReassignPoolImagewise(packs, edges_weight, reassign_opt);
      }

      int all_cluster_size = 0;
      for (size_t i = 0; i < packs.size(); ++i) {
        output_clusters[i] = packs[i].core_nodes;  // 最终聚合增量式的节点集合
        all_cluster_size += packs[i].core_nodes.size();
      }

      return true;
    }

    bool GeoCoviCluster::RunComputeOverlaps(
        const std::map<EdgeKey, int>& edges_weight,
        const OverlapOptions& overlap_opt,
        std::vector<std::vector<uint32_t>>& output_clusters) {

      auto overlaps =
          ComputePairwiseOverlaps(output_clusters, edges_weight, overlap_opt);
      for (size_t i = 0; i < output_clusters.size(); ++i) {
        for (size_t j = 0; j < output_clusters.size(); ++j) {
          if (i == j) continue;

          if (overlaps[i][j].size() != 0) {
            output_clusters[j].insert(output_clusters[j].end(),
                                      overlaps[i][j].begin(),
                                      overlaps[i][j].end());
          }
        }
      }

      return true;
    }

    bool GeoCoviCluster::RunAll(
        const std::map<EdgeKey, int>& covi_edges_weight,
        const std::map<EdgeKey, int>& geo_edges_weight,
        const uint32_t cluster_size,
        const uint32_t overlap_size,
        GeoCoviClusterOptions& opts,
        std::vector<std::vector<uint32_t>>& output_clusters) {
      std::map<EdgeKey, int> fusion_edges_weight;

      LOG(INFO) << "FuseGeoCoviEdgeWeight...";
      FuseGeoCoviEdgeWeight(covi_edges_weight,
                            geo_edges_weight,
                            opts.fusion_opts,
                            fusion_edges_weight);

      LOG(INFO) << "RunMetisSplit...";
      RunMetisSplit(
          fusion_edges_weight, cluster_size, 0, output_clusters);

      LOG(INFO) << "RunPruneCores...";
      std::vector<SubsetPack> packs;
      RunPruneCores(output_clusters, covi_edges_weight, opts.prune_opts, packs);

      LOG(INFO) << "RunReassignment...";
      RunReassignment(
          covi_edges_weight, opts.reassign_opts, packs, output_clusters);

      LOG(INFO) << "RunComputeOverlaps...";
      RunComputeOverlaps(covi_edges_weight, opts.overlap_opts, output_clusters);

      return true;
    }

    bool GeoCoviCluster::FillCluster(
        const std::vector<std::vector<uint32_t>>& output_clusters,
        colmap::SceneClustering::Cluster& root_cluster) {
      
        root_cluster.child_clusters.resize(output_clusters.size());

        for (int i = 0; i < output_clusters.size(); ++i) {
          colmap::SceneClustering::Cluster leaf_cluster;
          leaf_cluster.image_ids = output_clusters[i];
          root_cluster.child_clusters[i] = leaf_cluster;
        }

        // Sort child clusters by descending size of images and secondarily by
        // lowest
        // image id.
        std::sort(root_cluster.child_clusters.begin(),
                  root_cluster.child_clusters.end(),
                  [](const colmap::SceneClustering::Cluster& first,
                     const colmap::SceneClustering::Cluster& second) {
                    return first.image_ids.size() >= second.image_ids.size() &&
                           *std::min_element(first.image_ids.begin(),
                                             first.image_ids.end()) <
                               *std::min_element(second.image_ids.begin(),
                                                 second.image_ids.end());
                  });

        return true;
    }

    bool GeoCoviCluster::SavePLY_XYZ(const std::string& path,
                                     const std::vector<Eigen::Vector3d>& pts) {
      std::ofstream ofs(path);
      if (!ofs) return false;

      ofs << "ply\n";
      ofs << "format ascii 1.0\n";
      ofs << "element vertex " << pts.size() << "\n";
      ofs << "property float x\n";
      ofs << "property float y\n";
      ofs << "property float z\n";
      ofs << "end_header\n";

      ofs << std::setprecision(7);
      for (const auto& p : pts) {
        ofs << static_cast<float>(p.x()) << " " << static_cast<float>(p.y())
            << " " << static_cast<float>(p.z()) << "\n";
      }
      return true;
    }

    bool GeoCoviCluster::PrepareMetisData(
        const std::map<EdgeKey, int>& edges_weight,
        std::vector<EdgeKey>& vec_image_pairs,
        std::vector<int>& pair_weights) {
      vec_image_pairs.clear();
      pair_weights.clear();

      vec_image_pairs.reserve(edges_weight.size());
      pair_weights.reserve(edges_weight.size());
      for (const auto& kv : edges_weight) {
        vec_image_pairs.push_back(kv.first);
        pair_weights.push_back(static_cast<int>(kv.second));
      }

      return true;
    }

    bool GeoCoviCluster::FiltEdges(const std::map<EdgeKey, int>& input,
                                   std::map<EdgeKey, int>& output) {
      output.clear();

      auto pair_it = input.begin();
      while (pair_it != input.end()) {
        const auto& pair = pair_it->first;
        const auto& weight = pair_it->second;

        uint32_t img_i = pair.first;
        uint32_t img_j = pair.second;

        if (output.find({img_i, img_j}) == output.end() &&
            output.find({img_j, img_i}) == output.end()) {
          output[{img_i, img_j}] = weight;
          output[{img_j, img_i}] = weight;
        }

        ++pair_it;
      }

      return true;
    }

    bool GeoCoviCluster::MetisSplit(
        const std::vector<EdgeKey>& image_pairs,
        const std::vector<int> pair_weights,
        uint32_t cluster_size,
        uint32_t overlap_size,
        std::vector<std::vector<uint32_t>>& output_clusters) {
      output_clusters.clear();

      std::set<uint32_t> image_set;
      if (cluster_size <= 1) {
        for (auto& pair : image_pairs) {
          image_set.insert(pair.first);
          image_set.insert(pair.second);
        }

        std::vector<uint32_t> indice(image_set.begin(), image_set.end());
        output_clusters.push_back(std::move(indice));

        return true;
      }

      std::unique_ptr<colmap::SceneClustering::Cluster> root_cluster =
          std::make_unique<colmap::SceneClustering::Cluster>();

      root_cluster->image_ids.insert(
          root_cluster->image_ids.end(), image_set.begin(), image_set.end());

      colmap::SceneClustering::Options scene_cluster_options;
      scene_cluster_options.is_hierarchical = false;
      scene_cluster_options.image_overlap = 0;
      scene_cluster_options.branching = cluster_size;
      colmap::SceneClustering scene_clustering(scene_cluster_options);
      scene_clustering.Partition(image_pairs, pair_weights);
      
      auto leaf_clusters = scene_clustering.GetLeafClusters();

      // 整理输出结果
      for (int i = 0; i < leaf_clusters.size(); ++i) {
        std::vector<uint32_t> vec;
        for (auto image_id : leaf_clusters[i]->image_ids) {
          vec.push_back(image_id);
        }
        output_clusters.push_back(std::move(vec));
      }

      return true;
    }



    }  // namespace vgpart