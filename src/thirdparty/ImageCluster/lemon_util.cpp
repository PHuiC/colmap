#include "lemon_util.h"

#include <stack>
#include <functional>

#include "lemon/connectivity.h"

namespace vgpart {

SubsetPruneResult PruneSubsetByKcoreAndBridges(
    const std::vector<uint32_t>& subset_nodes,
    const std::map<EdgeKey, int>& edges_weight,
    int k_core,
    int min_edge_weight,
    bool detect_bridges) {
  SubsetPruneResult R;
  R.k_used = k_core;
  R.min_edge_weight = min_edge_weight;

  Graph g;
  lemon::ListGraph::NodeMap<uint32_t> node_id_map(g);

  std::unordered_map<uint32_t, Graph::Node> id2node;
  std::vector<Graph::Edge> edge_list;
  std::vector<EdgeKey> edge_keys;

  BuildGraphForSubset(subset_nodes,
                      edges_weight,
                      min_edge_weight,
                      g,
                      id2node,
                      node_id_map,
                      edge_list,
                      edge_keys);

  // k-core 剪枝
  KCorePrune(g, node_id_map, k_core, R.kept_nodes, R.removed_nodes);

  // 桥边诊断（对“保留子图”或原图都可以——这里对原图做）
  if (detect_bridges) {
    std::vector<EdgeKey> bridges;
    FindBridges(g, node_id_map, edge_list, edge_keys, bridges);
    R.bridges = std::move(bridges);
    R.has_bridge = !R.bridges.empty();
  }

  auto comps = ExtractConnectedComponents(
      R.kept_nodes, edges_weight, /*min_edge_weight=*/R.min_edge_weight);
  R.components = std::move(comps);
  if (!R.components.empty()) {
    size_t kept = R.kept_nodes.size();
    size_t lcc = R.components.front().size();
    R.lcc_ratio = kept ? (double)lcc / (double)kept : 1.0;
  }

  return R;
}
void BuildGraphForSubset(const std::vector<uint32_t>& subset_nodes,
                         const std::map<EdgeKey, int>& edges_weight,
                         int min_edge_weight,
                         Graph& g,
                         std::unordered_map<uint32_t, Graph::Node>& id2node,
                         lemon::ListGraph::NodeMap<uint32_t>& node_id_map,
                         std::vector<Graph::Edge>& edge_list,
                         std::vector<EdgeKey>& edge_keys) {
  id2node.clear();
  std::unordered_set<uint32_t> subset_set(subset_nodes.begin(),
                                          subset_nodes.end());

  // 加点
  for (uint32_t vid : subset_nodes) {
    auto n = g.addNode();
    id2node[vid] = n;
    node_id_map[n] = vid;
  }

  // 加边（仅子集内 & 权重达标）
  for (const auto& kv : edges_weight) {
    const auto& ek = kv.first;
    int w = kv.second;
    if (w < min_edge_weight) continue;
    if (subset_set.count(ek.first) && subset_set.count(ek.second)) {
      auto itA = id2node.find(ek.first);
      auto itB = id2node.find(ek.second);
      if (itA != id2node.end() && itB != id2node.end()) {
        auto e = g.addEdge(itA->second, itB->second);
        edge_list.push_back(e);
        edge_keys.push_back(ek);
      }
    }
  }
}

void KCorePrune(Graph& g,
                const lemon::ListGraph::NodeMap<uint32_t>& node_id_map,
                int k,
                std::vector<uint32_t>& kept,
                std::vector<uint32_t>& removed) {
  // 计算度
  lemon::ListGraph::NodeMap<int> deg(g, 0);
  std::vector<Graph::Node> nodes;
  for (Graph::NodeIt n(g); n != lemon::INVALID; ++n) {
    int d = 0;
    for (Graph::IncEdgeIt e(g, n); e != lemon::INVALID; ++e) ++d;
    deg[n] = d;
    nodes.push_back(n);
  }

  // 队列式剥皮
  std::stack<Graph::Node> st;
  lemon::ListGraph::NodeMap<bool> removed_flag(g, false);
  for (auto n : nodes)
    if (deg[n] < k) st.push(n);

  while (!st.empty()) {
    auto n = st.top();
    st.pop();
    if (removed_flag[n]) continue;
    removed_flag[n] = true;
    // 度更新
    for (Graph::IncEdgeIt e(g, n); e != lemon::INVALID; ++e) {
      auto u = g.oppositeNode(n, e);
      if (!removed_flag[u]) {
        deg[u] -= 1;
        if (deg[u] < k) st.push(u);
      }
    }
  }

  kept.clear();
  removed.clear();
  for (auto n : nodes) {
    uint32_t vid = node_id_map[n];
    if (removed_flag[n])
      removed.push_back(vid);
    else
      kept.push_back(vid);
  }
}

void FindBridges(const Graph& g,
                 const lemon::ListGraph::NodeMap<uint32_t>& node_id_map,
                 const std::vector<Graph::Edge>& edge_list,
                 const std::vector<EdgeKey>& edge_keys,
                 std::vector<EdgeKey>& bridges_out) {
  // 建索引：把 node 映射到 0..N-1，便于数组存取
  std::vector<Graph::Node> idx2node;
  idx2node.reserve(lemon::countNodes(g));
  lemon::ListGraph::NodeMap<int> node_idx(g, -1);
  int idx = 0;
  for (Graph::NodeIt n(g); n != lemon::INVALID; ++n) {
    node_idx[n] = idx++;
    idx2node.push_back(n);
  }
  const int N = idx;

  std::vector<int> tin(N, -1), low(N, -1), parent(N, -1);
  int timer = 0;
  std::vector<std::vector<std::pair<int, int>>> adj(N);  // (to, edge_index)

  // 边表
  for (int i = 0; i < (int)edge_list.size(); ++i) {
    auto e = edge_list[i];
    auto a = node_idx[g.u(e)];
    auto b = node_idx[g.v(e)];
    if (a < 0 || b < 0) continue;
    adj[a].push_back({b, i});
    adj[b].push_back({a, i});
  }

  std::vector<char> used_edge(edge_list.size(), 0);

  std::function<void(int, int)> dfs = [&](int v, int pe) {
    tin[v] = low[v] = timer++;
    for (auto [to, ei] : adj[v]) {
      if (ei == pe) continue;
      if (tin[to] != -1) {
        // 回边
        low[v] = std::min(low[v], tin[to]);
      } else {
        parent[to] = v;
        dfs(to, ei);
        low[v] = std::min(low[v], low[to]);
        if (low[to] > tin[v]) {
          // (v,to) 是桥
          used_edge[ei] = 1;
        }
      }
    }
  };

  for (int v = 0; v < N; ++v)
    if (tin[v] == -1) dfs(v, -1);

  bridges_out.clear();
  for (int i = 0; i < (int)edge_list.size(); ++i) {
    if (used_edge[i]) bridges_out.push_back(edge_keys[i]);
  }
}
std::vector<std::vector<uint32_t>> ExtractConnectedComponents(
    const std::vector<uint32_t>& kept_nodes,
    const std::map<EdgeKey, int>& edges_weight,
    int min_edge_weight) {
  Graph g;
  lemon::ListGraph::NodeMap<uint32_t> node_id_map(g);
  std::unordered_map<uint32_t, Graph::Node> id2node;
  id2node.reserve(kept_nodes.size());

  // 加点
  for (uint32_t vid : kept_nodes) {
    auto n = g.addNode();
    node_id_map[n] = vid;
    id2node[vid] = n;
  }

  // 加边（只在 kept_nodes 的诱导子图里，且权重大于等于阈值）
  std::unordered_set<uint32_t> kept_set(kept_nodes.begin(), kept_nodes.end());
  for (const auto& kv : edges_weight) {
    const auto& e = kv.first;
    int w = kv.second;
    if (w < min_edge_weight) continue;
    if (kept_set.count(e.first) && kept_set.count(e.second)) {
      g.addEdge(id2node[e.first], id2node[e.second]);
    }
  }

  // 连通分量
  lemon::ListGraph::NodeMap<int> comp(g, -1);
  int n_comp = lemon::connectedComponents(g, comp);

  // 收集各分量节点
  std::vector<std::vector<uint32_t>> comps(n_comp);
  for (Graph::NodeIt n(g); n != lemon::INVALID; ++n) {
    int cid = comp[n];
    comps[cid].push_back(node_id_map[n]);
  }

  // 按大小降序
  std::sort(comps.begin(), comps.end(), [](const auto& a, const auto& b) {
    return a.size() > b.size();
  });
  return comps;
}

long long CrossScoreComponentToSubset(
    const std::vector<uint32_t>& comp,
    const std::vector<uint32_t>& target_subset_nodes,
    const std::map<EdgeKey, int>& edges_weight,
    int min_edge_weight) {
  std::unordered_set<uint32_t> T(target_subset_nodes.begin(),
                                 target_subset_nodes.end());
  long long score = 0;
  for (uint32_t u : comp) {
    for (uint32_t v : target_subset_nodes) {
      auto it = edges_weight.find({u, v});
      if (it != edges_weight.end() && it->second >= min_edge_weight) {
        score += it->second;
      }
    }
  }
  return score;
}

bool KcoreFeasibleAfterMerge(const std::vector<uint32_t>& comp,
                             const std::vector<uint32_t>& target_subset_nodes,
                             const std::map<EdgeKey, int>& edges_weight,
                             int k_core,
                             int min_edge_weight,
                             double keep_ratio_thresh) {
  // 合并集合
  std::vector<uint32_t> union_nodes = target_subset_nodes;
  union_nodes.insert(union_nodes.end(), comp.begin(), comp.end());
  // 跑你已有的 k-core + 连接（可直接调用你封装的 prune 接口，或写一个 local
  // 版）
  auto R = PruneSubsetByKcoreAndBridges(union_nodes,
                                        edges_weight,
                                        k_core,
                                        min_edge_weight,
                                        /*detect_bridges=*/false);

  // 统计 comp 的保留比例
  std::unordered_set<uint32_t> kept(R.kept_nodes.begin(), R.kept_nodes.end());
  int kept_in_comp = 0;
  for (uint32_t u : comp)
    if (kept.count(u)) ++kept_in_comp;
  double ratio =
      comp.empty() ? 1.0 : (double)kept_in_comp / (double)comp.size();
  return ratio >= keep_ratio_thresh;
}

bool ReassignComponentsGreedy(std::vector<SubsetPack>& subsets,
                              const std::map<EdgeKey, int>& edges_weight,
                              const ReassignOptions& opt) {
  bool changed = false;

  // 对每个子集的每个 orphan 分量，找“最合适”的目标子集
  for (size_t s = 0; s < subsets.size(); ++s) {
    auto& src = subsets[s];
    std::vector<std::vector<uint32_t>> kept_orphans;  // 未能转移的，留下

    for (auto& comp : src.orphan_components) {
      long long best_score = std::numeric_limits<long long>::min();
      size_t best_t = SIZE_MAX;

      for (size_t t = 0; t < subsets.size(); ++t)
        if (t != s) {
          auto& dst = subsets[t];
          long long sc = CrossScoreComponentToSubset(
              comp, dst.core_nodes, edges_weight, opt.min_edge_weight);
          if (sc < opt.score_floor) continue;

          // 严格可行性检查
          if (!KcoreFeasibleAfterMerge(comp,
                                       dst.core_nodes,
                                       edges_weight,
                                       opt.k_core,
                                       opt.min_edge_weight,
                                       opt.keep_ratio_thresh))
            continue;

          if (sc > best_score) {
            best_score = sc;
            best_t = t;
          }
        }

      if (best_t != SIZE_MAX) {
        // 执行转移：合并到目标的 core 上（并重跑一次局部 k-core，保证稳）
        auto& dst = subsets[best_t];
        std::vector<uint32_t> union_nodes = dst.core_nodes;
        union_nodes.insert(union_nodes.end(), comp.begin(), comp.end());
        auto R = PruneSubsetByKcoreAndBridges(
            union_nodes, edges_weight, opt.k_core, opt.min_edge_weight, false);
        dst.core_nodes = std::move(
            R.kept_nodes);  // LCC之外的可再次归为
                            // orphan，但简单起见，这里只保留 k-core 结果
        changed = true;
      } else {
        kept_orphans.push_back(std::move(comp));  // 暂留
      }
    }

    src.orphan_components = std::move(kept_orphans);
  }

  return changed;
}

bool ReassignPoolImagewise(std::vector<SubsetPack>& subsets,
                           const std::map<EdgeKey, int>& edges_weight,
                           const ReassignOptions& opt) {
  bool changed = false;

  // 1) 收集所有 orphan 节点到池
  std::vector<uint32_t> pool;
  for (auto& s : subsets) {
    for (auto& comp : s.orphan_components) {
      pool.insert(pool.end(), comp.begin(), comp.end());
    }
    s.orphan_components.clear();
  }

  // 2) 逐图像挑目标子集（得分最大且可行）
  for (uint32_t u : pool) {
    long long best_score = std::numeric_limits<long long>::min();
    size_t best_t = SIZE_MAX;
    std::vector<uint32_t> single = {u};

    for (size_t t = 0; t < subsets.size(); ++t) {
      auto& dst = subsets[t];
      long long sc = CrossScoreComponentToSubset(
          single, dst.core_nodes, edges_weight, opt.min_edge_weight);

      if (sc <= opt.score_floor) continue;

      if (!KcoreFeasibleAfterMerge(single,
                                   dst.core_nodes,
                                   edges_weight,
                                   opt.k_core,
                                   opt.min_edge_weight,
                                   opt.keep_ratio_thresh))
        continue;

      if (sc > best_score) {
        best_score = sc;
        best_t = t;
      }
    }

    if (best_t != SIZE_MAX) {
      auto& dst = subsets[best_t];
      std::vector<uint32_t> union_nodes = dst.core_nodes;
      union_nodes.push_back(u);
      auto R = PruneSubsetByKcoreAndBridges(
          union_nodes, edges_weight, opt.k_core, opt.min_edge_weight, false);
      dst.core_nodes = std::move(R.kept_nodes);
      changed = true;
    } else {
      // 3) 实在没地方可放，可选择：新建子集 / 丢弃 / 暂存回 orphan
      // 这里给个保守策略：回到“最相近”的子集作为 orphan，便于下一轮组合尝试
      size_t t_fallback = SIZE_MAX;
      long long best_any = std::numeric_limits<long long>::min();
      for (size_t t = 0; t < subsets.size(); ++t) {
        long long sc = CrossScoreComponentToSubset(
            single, subsets[t].core_nodes, edges_weight, opt.min_edge_weight);
        if (sc > best_any) {
          best_any = sc;
          t_fallback = t;
        }
      }
      if (t_fallback != SIZE_MAX) {
        subsets[t_fallback].orphan_components.push_back({u});
      }
    }
  }

  return changed;
}  

int EffectiveDegreeToSet(uint32_t u,
                         const std::unordered_set<uint32_t>& T,
                         const std::map<EdgeKey, int>& edges_weight,
                         int min_edge_weight) {
  int deg = 0;
  for (uint32_t v : T) {
    auto it_uv = edges_weight.find({u, v});
    if (it_uv != edges_weight.end() && it_uv->second >= min_edge_weight) ++deg;

    auto it_vu = edges_weight.find({v, u});
    if (it_vu != edges_weight.end() && it_vu->second >= min_edge_weight) ++deg;
  }
  return deg;
}

long long CrossScoreToSet(uint32_t u,
                          const std::unordered_set<uint32_t>& T,
                          const std::map<EdgeKey, int>& edges_weight,
                          int min_edge_weight) {
  long long s = 0;
  for (uint32_t v : T) {
    auto it_uv = edges_weight.find({u, v});
    if (it_uv != edges_weight.end() && it_uv->second >= min_edge_weight)
      s += it_uv->second;

    auto it_vu = edges_weight.find({v, u});
    if (it_vu != edges_weight.end() && it_vu->second >= min_edge_weight)
      s += it_vu->second;
  }
  return s;
}

std::vector<uint32_t> SelectOverlapFromSrcToDst(
    const std::vector<uint32_t>& A_core,
    const std::vector<uint32_t>& B_core,
    const std::map<EdgeKey, int>& edges_weight,
    const OverlapOptions& p) {
  // 目标集合：从 B_core 开始，逐步加入选中的 A 的节点（增量判定）
  std::unordered_set<uint32_t> target(B_core.begin(), B_core.end());

  // 候选：A_core 去掉与 B_core 的交集（避免把 B 已有的点再选一次）
  std::unordered_set<uint32_t> Bset(B_core.begin(), B_core.end());
  std::vector<uint32_t> candidates;
  candidates.reserve(A_core.size());
  for (uint32_t u : A_core)
    if (!Bset.count(u)) candidates.push_back(u);

  // 相关性打分：先按（对 B_core 的）跨边权重和排序，从大到小
  std::vector<std::pair<uint32_t, long long>> scored;
  scored.reserve(candidates.size());
  for (uint32_t u : candidates) {
    long long sc = CrossScoreToSet(
        u,
        std::unordered_set<uint32_t>(B_core.begin(), B_core.end()),
        edges_weight,
        p.min_edge_weight);
    scored.emplace_back(u, sc);
  }
  std::sort(scored.begin(), scored.end(), [](const auto& a, const auto& b) {
    return (a.second != b.second) ? a.second > b.second : a.first < b.first;
  });

  // 贪心选择：必须在“当前 target 集合”上度数 >= k_core
  std::vector<uint32_t> picked;
  picked.reserve(std::min<int>((int)scored.size(), p.K));
  for (const auto& pr : scored) {
    if ((int)picked.size() >= p.K) break;
    uint32_t u = pr.first;
    // 度判定：注意度是针对“当前 target（B_core ∪ 已选 A）”
    int deg = EffectiveDegreeToSet(u, target, edges_weight, p.min_edge_weight);
    if (deg >= p.k_core) {
      picked.push_back(u);
      target.insert(u);  // 加入增强集合，有助于后续候选过线
    }
  }
  return picked;
}

std::vector<std::vector<std::vector<uint32_t>>> ComputePairwiseOverlaps(
    const std::vector<std::vector<uint32_t>>& cores,
    const std::map<EdgeKey, int>& edges_weight,
    const OverlapOptions& p) {
  const size_t N = cores.size();
  std::vector<std::vector<std::vector<uint32_t>>> overlaps(
      N, std::vector<std::vector<uint32_t>>(N));

  // 1) 计算每个 i 与其它子集的相关度，并选前 top_blocks 作为“种子方向”
  std::vector<std::vector<size_t>> topJs(N);
  for (size_t i = 0; i < N; ++i) {
    std::vector<std::pair<long long, size_t>> scored;
    scored.reserve(N - 1);
    for (size_t j = 0; j < N; ++j) {
      if (i == j) continue;
      long long sc = CrossScoreComponentToSubset(
          cores[i], cores[j], edges_weight, p.min_edge_weight);
      scored.emplace_back(sc, j);
    }
    std::sort(scored.begin(), scored.end(), [](const auto& a, const auto& b) {
      if (a.first != b.first) return a.first > b.first;
      return a.second < b.second;
    });

    int take = (p.top_blocks <= 0)
                   ? (int)scored.size()
                   : std::min<int>(p.top_blocks, (int)scored.size());
    topJs[i].reserve(take);
    for (int k = 0; k < take; ++k) {
      // 相关度>0 才作为“种子”，避免无意义计算
      if (scored[k].first <= 0) break;
      topJs[i].push_back(scored[k].second);
    }
  }

  // 2) 对称性保护：做对称闭包，不受 top-n 限制
  std::vector<std::vector<char>> sel(N, std::vector<char>(N, 0));
  for (size_t i = 0; i < N; ++i) {
    for (size_t j : topJs[i]) {
      sel[i][j] = 1;  // 种子方向
      sel[j][i] = 1;  // ★ 对称闭包（即使 j 的 top-n 里没有 i 也强制开启）
    }
  }

  for (size_t i = 0; i < N; ++i) {
    for (size_t j = 0; j < N; ++j) {
      if (i == j || !sel[i][j]) continue;

      std::vector<uint32_t> A_core = cores[i];
      std::vector<uint32_t> B_core = cores[j];
      std::vector<uint32_t> overlap_result;
      for (int i = 0; i < 10; ++i) {
        A_core.insert(
            A_core.end(), overlap_result.begin(), overlap_result.end());
        auto overlaps =
            SelectOverlapFromSrcToDst(B_core, A_core, edges_weight, p);
        overlap_result.insert(
            overlap_result.end(), overlaps.begin(), overlaps.end());
        if (overlap_result.size() >= p.K) {
          overlap_result.resize(p.K);
          break;
        }
      }
      overlaps[j][i] = overlap_result;
    }
  }
  return overlaps;
}

DegreeResult ScoreViewGraphDegrees(
    const std::map<EdgeKey, int>& covi_edges_weight, int min_weight) {

    DegreeResult result;

    // 1) LEMON 建图
    lemon::ListGraph G;
    lemon::ListGraph::NodeMap<uint32_t> id_of_node(G);  // Node -> node_id
    std::unordered_map<uint32_t, lemon::ListGraph::Node> node_of_id;

    auto ensure_node = [&](uint32_t id) -> lemon::ListGraph::Node {
      auto it = node_of_id.find(id);
      if (it != node_of_id.end()) return it->second;
      auto n = G.addNode();
      node_of_id.emplace(id, n);
      id_of_node[n] = id;
      return n;
    };

    for (const auto& kv : covi_edges_weight) {
      const uint32_t u = kv.first.first;
      const uint32_t v = kv.first.second;
      if (u == v) continue;
      auto nu = ensure_node(u);
      auto nv = ensure_node(v);
      G.addEdge(nu, nv);
    }

    // 过滤后若没有任何节点，按约定返回 0/false
    if (lemon::countNodes(G) == 0) {
      result.stats = DegreeStats{0, 0, 0.0, 0.0, false};
      result.node_degree.clear();
      return result;
    }

    // 2) 统计每个节点的度（LEMON: countIncEdges）
    uint32_t min_deg = std::numeric_limits<uint32_t>::max();
    uint32_t max_deg = 0;
    uint64_t sum_deg = 0;

    std::vector<uint32_t> deg_vals;
    deg_vals.reserve(lemon::countNodes(G));

    for (lemon::ListGraph::NodeIt n(G); n != lemon::INVALID; ++n) {
      const uint32_t d = static_cast<uint32_t>(lemon::countIncEdges(G, n));
      const uint32_t id = id_of_node[n];
      result.node_degree[id] = d;

      min_deg = std::min(min_deg, d);
      max_deg = std::max(max_deg, d);
      sum_deg += d;
      deg_vals.push_back(d);
    }
    std::sort(deg_vals.begin(), deg_vals.end());

    const double avg =
        static_cast<double>(sum_deg) / static_cast<double>(deg_vals.size());
    const double med = MedianOfSorted(deg_vals);

    // 3) 连通性（LEMON: connectedComponents）
    lemon::ListGraph::NodeMap<int> comp(G);
    const int num_comp = lemon::connectedComponents(G, comp);
    const bool single = (num_comp == 1);

    result.stats.min_degree = min_deg;
    result.stats.max_degree = max_deg;
    result.stats.avg_degree = avg;
    result.stats.median_degree = med;
    result.stats.is_single_component = single;
    return result;

}
}  // namespace vgpart
