// triaccel.cpp
//
// このファイルは triaccel の C++ バックエンド実装（pybind11 バインディング）です。
// 役割：
//  - サイト条件に基づくイベント（天球上の方向ベクトル）の乱数生成
//  - 角距離しきい値 d_deg 以内の隣接関係の構築
//  - k-クリーク（完全グラフ）のカウント／必要に応じて三重項のヒスト更新
//  - （debug 時）テキストログの出力と段階ログの標準エラー出力
//  - Python dict で counts（およびヒスト／エッジ配列）を返却
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <vector>
#include <cmath>
#include <cstdint>
#include <random>
#include <chrono>
#include <thread>
#include <mutex>
#include <atomic>
#include <cstdio> // fprintf, fputc
#include <string>
#include <fstream>
#include <functional>

#include "detail/utils.hpp"
// グラフ処理（隣接構築・クリークカウント）
#include "detail/graph_ops.hpp"

namespace py = pybind11;



// 事前計算済みの RA/Dec 配列を使ってイベント・ヒストグラムを更新（再計算を回避）
static void update_event_histograms_if_needed(
    int N,
    const std::vector<double> &ra_deg_arr,
    const std::vector<double> &dec_deg_arr,
    int tid,
    bool return_histograms,
    int bins_ra,
    const std::function<void(double, double, int &, int &)> &bin2d,
    std::vector<std::vector<long long>> &Evt2d_loc)
{
    if (!return_histograms)
    {
        return;
    }

    for (int i = 0; i < N; ++i)
    {
        int ir, id;
        bin2d(ra_deg_arr[i], dec_deg_arr[i], ir, id);
        Evt2d_loc[tid][(size_t)id * (size_t)bins_ra + (size_t)ir] += 1;
    }
}

// --- end helpers ---

// 乱数生成器（[0,1), [0,2π) を手軽に取得）
struct RNG
{
    std::mt19937_64 eng;
    std::uniform_real_distribution<double> U01, U02pi;
    RNG(uint64_t seed) : eng(seed), U01(0.0, 1.0), U02pi(0.0, 2 * PI) {}
    inline double u01() { return U01(eng); }
    inline double u02pi() { return U02pi(eng); }
};

// サイト（観測点）パラメタの内部表現（前計算済み三角関数を保持）
struct SiteC
{
    double sphi, cphi, slam, clam, smax;
    int n;
};
// サイト1件の生成（緯度・経度・最大天頂角 zmax と、イベント数 n）
static SiteC make_site(double lat_deg, double lon_deg, double zmax_deg, int n)
{
    double phi = rad(lat_deg), lam = rad(lon_deg);
    return SiteC{std::sin(phi), std::cos(phi), std::sin(lam), std::cos(lam),
                 std::sin(rad(zmax_deg)), n};
}
// Python の sites 配列（辞書の列）を SiteC ベクタに変換
static std::vector<SiteC> parse_sites(py::sequence sites)
{
    std::vector<SiteC> v;
    v.reserve(py::len(sites));
    for (auto item : sites)
    {
        auto d = py::cast<py::dict>(item);
        int n = d["n"].cast<int>();
        if (n <= 0)
            continue;
        v.push_back(make_site(d["lat"].cast<double>(),
                              d["lon"].cast<double>(),
                              d["zmax"].cast<double>(), n));
    }
    return v;
}
// 全サイトのイベント総数 N を返す
static int total_events(const std::vector<SiteC> &S)
{
    int N = 0;
    for (auto &s : S)
        N += s.n;
    return N;
}

// 1サイトぶんのイベント方向を ECI 座標で生成
// ENU → ECEF → ECI（地球自転の Z 回転）で変換
static void gen_site_eci(const SiteC &S, int offset, RNG &rng,
                         std::vector<double> &x, std::vector<double> &y, std::vector<double> &z)
{
    for (int i = 0; i < S.n; ++i)
    {
        double lst = rng.u02pi();
        double clst = std::cos(lst), slst = std::sin(lst);
        double az = rng.u02pi();
        double caz = std::cos(az), saz = std::sin(az);
        // θ ~ sinθ cosθ : sinθ=smax*sqrt(u), cosθ=sqrt(1-smax^2 u)
        double u = rng.u01();
        double sinth = S.smax * std::sqrt(u);
        double costh = std::sqrt(std::max(0.0, 1.0 - (S.smax * S.smax) * u));
        // ENU（alt=π/2-θ → cos(alt)=sinθ, sin(alt)=cosθ）
        // 水平成分 = cos(alt) = sinθ、上向き成分 = sin(alt) = cosθ
        double e = sinth * saz;
        double nN = sinth * caz;
        double uU = costh;
        // ENU→ECEF
        double xe = e * (-S.slam) + nN * (-S.sphi * S.clam) + uU * (S.cphi * S.clam);
        double ye = e * (S.clam) + nN * (-S.sphi * S.slam) + uU * (S.cphi * S.slam);
        double ze = nN * (S.cphi) + uU * (S.sphi);
        // ECEF→ECI（Z回転=LST）
        int k = offset + i;
        x[k] = xe * clst - ye * slst;
        y[k] = xe * slst + ye * clst;
        z[k] = ze;
    }
}

// N<=128向け：2×64bit のビットマスクで隣接行列を構築
// cos_thr 以上の内積（=角距離がしきい値以内）なら辺を立てる
// ---- helpers ----

// シミュレーション本体
// 引数：
//  - M: 試行回数
//  - sites: サイト辞書の列（lat/lon/zmax/n）
//  - d_deg: 角距離しきい値（度）
//  - seed/bins/threads/progress: 各種オプション
//  - debug: デバッグ（ログ出力・Mクランプ・段階ログ）
//  - cluster_size: k-クリークの k（2以上）

// --- ワーカー用バッファ（各試行で使い回すメモリ領域） ---
struct ThreadBuf
{
    std::vector<double>   x, y, z;
    std::vector<uint64_t> m0, m1;
    std::vector<int>      site_id;
    std::vector<char>     in_cluster;
    std::vector<double>   ra_deg_arr, dec_deg_arr;
    ThreadBuf(int N)
        : x(N), y(N), z(N), m0(N), m1(N), site_id(N, 0), in_cluster(N, 0),
          ra_deg_arr(N), dec_deg_arr(N) {}
    void clear_cluster() { std::fill(in_cluster.begin(), in_cluster.end(), 0); }
};

// 乱数生成によりイベント方向（ECI）を生成し、RA/Dec を配列にキャッシュする
static void generate_events_and_cache(
    const std::vector<SiteC> &S,
    RNG &rng,
    ThreadBuf &B)
{
    int off = 0;
    for (size_t sid = 0; sid < S.size(); ++sid)
    {
        const auto &site = S[sid];
        gen_site_eci(site, off, rng, B.x, B.y, B.z);
        for (int i = 0; i < site.n; ++i)
        {
            B.site_id[off + i] = (int)sid;
        }
        off += site.n;
    }
    const int N = (int)B.x.size();
    for (int i = 0; i < N; ++i)
    {
        double ra_deg, dec_deg;
        compute_ra_dec_deg_from_xyz(B.x[i], B.y[i], B.z[i], ra_deg, dec_deg);
        B.ra_deg_arr[i] = ra_deg;
        B.dec_deg_arr[i] = dec_deg;
        B.in_cluster[i] = 0;
    }
}

// debug 用：N<=128 の k-クリーク列挙をテキストに書き出す（i<j<... を保証）
static void dump_k_cliques_u128(
    int K,
    const std::vector<uint64_t> &m0,
    const std::vector<uint64_t> &m1,
    int N,
    std::ofstream &cf)
{
    if (K < 2 || N <= 0) return;

    auto and_with_nb = [&](uint64_t c0, uint64_t c1, int v, uint64_t &o0, uint64_t &o1) {
        o0 = c0 & m0[v];
        o1 = c1 & m1[v];
    };
    auto mask_gt = [&](int v, uint64_t &c0, uint64_t &c1) {
        // インデックス > v だけ残す
        if (v < 63) {
            uint64_t mask0 = ~((1ULL << (v + 1)) - 1ULL);
            c0 &= mask0;
            // c1 はそのまま
        } else if (v == 63) {
            c0 = 0ULL; // 0..63 は全て落とす
        } else { // v>=64
            c0 = 0ULL;
            unsigned vb = (unsigned)(v - 64);
            if (vb < 63) {
                uint64_t mask1 = ~((1ULL << (vb + 1)) - 1ULL);
                c1 &= mask1;
            } else {
                c1 = 0ULL;
            }
        }
    };

    std::vector<int> cur;
    cur.reserve(K);

    std::function<void(int,int,uint64_t,uint64_t)> dfs = [&](int last, int depth, uint64_t c0, uint64_t c1) {
        if (depth == K) {
            // 出力
            for (int i = 0; i < K; ++i) {
                if (i) cf << '\t';
                cf << cur[i];
            }
            cf << '\n';
            return;
        }
        // 枝刈り（候補の個数 < 残り必要数）
        int remain = K - depth;
        int cand_cnt = (int)pop64(c0) + (int)pop64(c1);
        if (cand_cnt < remain) return;
        // 低位ビットから順に展開
        uint64_t t0 = c0, t1 = c1;
        while (t0) {
            unsigned b = (unsigned)__builtin_ctzll(t0);
            int v = (int)b;
            uint64_t n0, n1; and_with_nb(c0, c1, v, n0, n1);
            mask_gt(v, n0, n1);
            cur.push_back(v);
            dfs(v, depth + 1, n0, n1);
            cur.pop_back();
            t0 &= (t0 - 1ULL);
        }
        while (t1) {
            unsigned b = (unsigned)__builtin_ctzll(t1);
            int v = 64 + (int)b;
            if (v >= N) break;
            uint64_t n0, n1; and_with_nb(c0, c1, v, n0, n1);
            mask_gt(v, n0, n1);
            cur.push_back(v);
            dfs(v, depth + 1, n0, n1);
            cur.pop_back();
            t1 &= (t1 - 1ULL);
        }
    };

    // 先頭頂点を 0..N-1 で回す（候補は近傍に限定）
    for (int v = 0; v < N; ++v) {
        uint64_t c0 = m0[v];
        uint64_t c1 = m1[v];
        mask_gt(v, c0, c1);
        cur.clear(); cur.push_back(v);
        dfs(v, 1, c0, c1);
    }
}

// N<=128 の高速パス（ビットマスク 2 語）で 1 試行を処理
static void process_smallN_trial(
    int t,
    int M,
    int N,
    int tid,
    int cluster_size,
    double cos_thr,
    int bins_ra,
    bool debug,
    bool return_histograms,
    const std::function<void(double, double, int &, int &)> &bin2d,
    const std::string &log_text_dir,
    ThreadBuf &B,
    std::vector<std::vector<long long>> &Evt2d_loc,
    std::vector<std::vector<long long>> &Tri2d_loc,
    std::vector<int32_t> &counts_vec)
{
    build_masks_u128(B.x, B.y, B.z, cos_thr, B.m0, B.m1);

    if (debug)
    {
        long long edges = 0;
        int dmin = N;
        int dmax = 0;
        double dsum = 0.0;
        for (int i = 0; i < N; ++i)
        {
            int di = pop64(B.m0[i]) + pop64(B.m1[i]);
            dsum += di;
            if (di < dmin) dmin = di;
            if (di > dmax) dmax = di;
        }
        edges = (long long)(dsum / 2.0);
        std::fprintf(stderr,
                     "t=%d/%d build adjacency (edges=%lld, deg[min/avg/max]=%d/%.2f/%d)\n",
                     t + 1, M, edges, dmin, (N ? dsum / N : 0.0), dmax);
    }

    if (cluster_size == 3)
    {
        long long tri_cnt = count_triangles_u128(B.m0, B.m1);
        counts_vec[t] = (int32_t)tri_cnt;

        if (return_histograms)
        {
            enumerate_triangles_hist_from_masks_u128(
                N, tid,
                B.x, B.y, B.z,
                B.m0, B.m1,
                /*mark_members=*/debug,
                B.in_cluster,
                /*return_hist2d=*/return_histograms,
                bins_ra,
                bin2d,
                Tri2d_loc,
                /*tri2d_push_centroid=*/true);
        }
    }
    else if (cluster_size == 2)
    {
        long long pair_cnt = count_pairs_u128(B.m0, B.m1, N, &B.in_cluster, debug);
        counts_vec[t] = (int32_t)pair_cnt;
        if (return_histograms)
        {
            enumerate_pairs_and_push_hist_u128(N, tid, B.m0, B.m1,
                                               B.x, B.y, B.z,
                                               /*return_hist2d=*/return_histograms,
                                               bins_ra,
                                               bin2d,
                                               Tri2d_loc);
        }
    }
    else
    {
        long long kcnt = count_k_cliques_u128(cluster_size, B.m0, B.m1, N, &B.in_cluster, debug);
        counts_vec[t] = (int32_t)kcnt;
        // For k>=4, also push centroid histogram if needed
        if (return_histograms)
        {
            int W_nb = 0;
            std::vector<uint64_t> NB;
            masks_to_NB(B.m0, B.m1, N, NB, W_nb);
            long long dummy = count_k_cliques_bitadj(
                cluster_size, NB, N, W_nb, &B.in_cluster, debug,
                tid, &B.x, &B.y, &B.z,
                /*return_hist2d=*/return_histograms,
                bins_ra,
                &bin2d,
                &Tri2d_loc);
            (void)dummy;
        }
    }


    if (debug)
    {
        std::fprintf(stderr,
                     "t=%d/%d count k-cliques (k=%d, fastpath=%s)\n",
                     t + 1, M, cluster_size, (cluster_size == 3 ? "true" : "false"));
        try
        {
            char evname[64];
            std::snprintf(evname, sizeof(evname), "sim_%05d_events.txt", t + 1);
            char clqname[64];
            std::snprintf(clqname, sizeof(clqname), "sim_%05d_cliques_k%d.txt", t + 1, cluster_size);
            std::string evpath = log_text_dir + "/" + std::string(evname);
            std::string clqpath = log_text_dir + "/" + std::string(clqname);
            std::ofstream evf(evpath.c_str());
            evf << "# idx\tsite\tra_deg\tdec_deg\tin_cluster\n";
            for (int i = 0; i < N; ++i)
            {
                double ra_wr = wrap_deg180(B.ra_deg_arr[i]);
                evf << i << '\t' << B.site_id[i] << '\t' << ra_wr
                    << '\t' << B.dec_deg_arr[i] << '\t' << (int)B.in_cluster[i] << '\n';
            }
            evf.close();
            std::ofstream cf(clqpath.c_str());
            cf << "# v0\tv1\tv2...\n";
            if (cluster_size == 3)
            {
                auto has_edge = [&](int a, int b) -> bool {
                    if (b < 64) return (B.m0[a] >> b) & 1ULL;
                    else        return (B.m1[a] >> (b - 64)) & 1ULL;
                };
                for (int i0 = 0; i0 < N; ++i0)
                {
                    for (int j0 = i0 + 1; j0 < N; ++j0)
                    {
                        if (!has_edge(i0, j0)) continue;
                        uint64_t c0 = B.m0[i0] & B.m0[j0];
                        uint64_t c1 = B.m1[i0] & B.m1[j0];
                        if (j0 < 64)
                        {
                            uint64_t mask_lo = (j0 < 63) ? ~((1ULL << (j0 + 1)) - 1ULL) : 0ULL;
                            c0 &= mask_lo;
                        }
                        else
                        {
                            c0 = 0ULL;
                            unsigned jb = (unsigned)(j0 - 64);
                            uint64_t mask_hi = (jb < 63) ? ~((1ULL << (jb + 1)) - 1ULL) : 0ULL;
                            c1 &= mask_hi;
                        }
                        while (c0)
                        {
                            unsigned b = (unsigned)__builtin_ctzll(c0);
                            int k = (int)b;
                            if (k < N) { cf << i0 << '\t' << j0 << '\t' << k << '\n'; }
                            c0 &= (c0 - 1ULL);
                        }
                        while (c1)
                        {
                            unsigned b = (unsigned)__builtin_ctzll(c1);
                            int k = 64 + (int)b;
                            if (k < N) { cf << i0 << '\t' << j0 << '\t' << k << '\n'; }
                            c1 &= (c1 - 1ULL);
                        }
                    }
                }
            }
            else if (cluster_size == 2)
            {
                auto has_edge = [&](int a, int b) -> bool {
                    if (b < 64) return (B.m0[a] >> b) & 1ULL;
                    else        return (B.m1[a] >> (b - 64)) & 1ULL;
                };
                // list all edges (i<j)
                for (int i0 = 0; i0 < N; ++i0)
                {
                    for (int j0 = i0 + 1; j0 < N; ++j0)
                    {
                        if (has_edge(i0, j0)) { cf << i0 << '\t' << j0 << '\n'; }
                    }
                }
            }
            else
            {
                // k>=4 の列挙（u128 マスク）
                cf << "# v0\tv1\t... (k=" << cluster_size << ")\n";
                dump_k_cliques_u128(cluster_size, B.m0, B.m1, N, cf);
            }
            cf.close();
        }
        catch (...) { }
    }
}

// debug 用：一般幅 NB から k-クリーク列挙をテキストに書き出す（i<j<... を保証）
static void dump_k_cliques_bitadj(
    int K,
    const std::vector<uint64_t> &NB,
    int N,
    int W,
    std::ofstream &cf)
{
    if (K < 2 || N <= 0) return;

    auto and_row = [&](const std::vector<uint64_t> &a, int v, std::vector<uint64_t> &o) {
        const uint64_t *row = &NB[(size_t)v * (size_t)W];
        for (int w = 0; w < W; ++w) o[w] = a[w] & row[w];
    };
    auto mask_gt = [&](int v, std::vector<uint64_t> &c) {
        int wb = v >> 6; int bb = v & 63;
        for (int w = 0; w < W; ++w) {
            if (w < wb) c[w] = 0ULL;
            else if (w == wb) {
                uint64_t mask = (bb < 63) ? ~((1ULL << (bb + 1)) - 1ULL) : 0ULL;
                c[w] &= mask;
            }
        }
        // 上位語はそのまま
    };
    auto popcnt_vec = [&](const std::vector<uint64_t> &a) {
        int s = 0; for (int w = 0; w < W; ++w) s += (int)pop64(a[w]); return s;
    };

    std::vector<int> cur; cur.reserve(K);
    std::vector<uint64_t> cand(W), next(W);

    std::function<void(int,int,const std::vector<uint64_t>&)> dfs =
    [&](int last, int depth, const std::vector<uint64_t> &c) {
        if (depth == K) {
            for (int i = 0; i < K; ++i) { if (i) cf << '\t'; cf << cur[i]; } cf << '\n';
            return;
        }
        int remain = K - depth; if (popcnt_vec(c) < remain) return;
        // 語ごとに走査
        for (int w = 0; w < W; ++w) {
            uint64_t bits = c[w];
            while (bits) {
                unsigned b = (unsigned)__builtin_ctzll(bits);
                int v = (w << 6) + (int)b; if (v >= N) break;
                // 次の候補 = c ∧ NB[v]、かつ >v に制限
                and_row(c, v, next);
                mask_gt(v, next);
                cur.push_back(v);
                dfs(v, depth + 1, next);
                cur.pop_back();
                bits &= (bits - 1ULL);
            }
        }
    };

    // 先頭頂点 v をループ（候補 = NB[v] のみ）
    for (int v = 0; v < N; ++v) {
        for (int w = 0; w < W; ++w) cand[w] = NB[(size_t)v * (size_t)W + (size_t)w];
        mask_gt(v, cand);
        cur.clear(); cur.push_back(v);
        dfs(v, 1, cand);
    }
}

// N>128 の一般パス（任意幅の語配列）で 1 試行を処理
static void process_largeN_trial(
    int t,
    int M,
    int N,
    int tid,
    int cluster_size,
    double cos_thr,
    int bins_ra,
    bool debug,
    bool return_histograms,
    const std::function<void(double, double, int &, int &)> &bin2d,
    const std::string &log_text_dir,
    ThreadBuf &B,
    std::vector<std::vector<long long>> &Evt2d_loc,
    std::vector<std::vector<long long>> &Tri2d_loc,
    std::vector<int32_t> &counts_vec)
{
    int W = (N + 63) / 64;
    std::vector<uint64_t> NB((size_t)N * (size_t)W, 0ULL);
    auto set_edge = [&](int a, int b)
    {
        int wb = b >> 6; int bb = b & 63;
        NB[(size_t)a * (size_t)W + (size_t)wb] |= (1ULL << bb);
    };
    // ビット隣接 NB の set/test はローカルラムダで十分に明瞭

    long long edges_undirected = 0;
    std::vector<int> degv(N, 0);
    for (int i2 = 0; i2 < N; ++i2)
    {
        double xi = B.x[i2], yi = B.y[i2], zi = B.z[i2];
        for (int j2 = i2 + 1; j2 < N; ++j2)
        {
            double dot = xi * B.x[j2] + yi * B.y[j2] + zi * B.z[j2];
            if (dot >= cos_thr)
            {
                set_edge(i2, j2); set_edge(j2, i2);
                edges_undirected++; degv[i2]++; degv[j2]++;
            }
        }
    }

    if (debug)
    {
        double sum = 0.0; int dmin = N, dmax = 0;
        for (int di : degv) { sum += di; if (di < dmin) dmin = di; if (di > dmax) dmax = di; }
        std::fprintf(stderr,
                     "t=%d/%d build adjacency (edges=%lld, deg[min/avg/max]=%d/%.2f/%d)\n",
                     t + 1, M, edges_undirected, dmin, (N ? sum / N : 0.0), dmax);
    }

    if (cluster_size == 3)
    {
        long long tri_cnt_gen = count_and_hist_triangles_bitadj(
            NB, N, W, tid,
            B.x, B.y, B.z,
            debug, B.in_cluster,
            return_histograms,
            bins_ra,
            bin2d,
            Tri2d_loc,
            /*tri2d_push_centroid=*/true);
        counts_vec[t] = (int32_t)tri_cnt_gen;
    }
    else if (cluster_size == 2)
    {
        // We already computed the undirected edge count while building NB.
        // Use it directly for k=2 to keep smallN/largeN paths consistent.
        long long pair_cnt = edges_undirected;
        if (debug)
        {
            long long check = count_pairs_bitadj(NB, N, W, &B.in_cluster, debug);
            if (check != pair_cnt)
            {
                std::fprintf(stderr, "[warn] pair_cnt mismatch: edges_undirected=%lld, recount=%lld\n",
                             pair_cnt, check);
            }
        }
        counts_vec[t] = (int32_t)pair_cnt;
        if (return_histograms)
        {
            enumerate_pairs_and_push_hist_bitadj(NB, N, W, tid,
                                                 B.x, B.y, B.z,
                                                 /*return_hist2d=*/return_histograms,
                                                 bins_ra,
                                                 bin2d,
                                                 Tri2d_loc);
        }
    }
    else
    {
        long long kcnt = count_k_cliques_bitadj(
            cluster_size, NB, N, W, &B.in_cluster, debug,
            tid, &B.x, &B.y, &B.z,
            return_histograms,
            bins_ra,
            &bin2d,
            &Tri2d_loc);
        counts_vec[t] = (int32_t)kcnt;
    }

    if (debug)
    {
        std::fprintf(stderr,
                     "t=%d/%d count k-cliques (k=%d, fastpath=%s)\n",
                     t + 1, M, cluster_size, (cluster_size == 3 ? "true" : "false"));
        try
        {
            char evname[64]; std::snprintf(evname, sizeof(evname), "sim_%05d_events.txt", t + 1);
            std::string evpath = log_text_dir + "/" + std::string(evname);
            std::ofstream evf(evpath.c_str());
            evf << "# idx\tsite\tra_deg\tdec_deg\tin_cluster\n";
            for (int i = 0; i < N; ++i)
            {
                double ra_wr = wrap_deg180(B.ra_deg_arr[i]);
                evf << i << '\t' << B.site_id[i] << '\t' << ra_wr
                    << '\t' << B.dec_deg_arr[i] << '\t' << (int)B.in_cluster[i] << '\n';
            }
            evf.close();
            // also dump cliques in debug mode
            char clqname[64]; std::snprintf(clqname, sizeof(clqname), "sim_%05d_cliques_k%d.txt", t + 1, cluster_size);
            std::string clqpath = log_text_dir + "/" + std::string(clqname);
            std::ofstream cf(clqpath.c_str());
            if (cluster_size == 3)
            {
                cf << "# v0\tv1\tv2...\n";
                for (int i0 = 0; i0 < N; ++i0)
                {
                    for (int j0 = i0 + 1; j0 < N; ++j0)
                    {
                        int wb = j0 >> 6; int bb = j0 & 63;
                        if ((NB[(size_t)i0 * (size_t)W + (size_t)wb] & (1ULL << bb)) == 0ULL) continue;
                        for (int w = wb; w < W; ++w)
                        {
                            uint64_t c = NB[(size_t)i0 * (size_t)W + (size_t)w] & NB[(size_t)j0 * (size_t)W + (size_t)w];
                            if (w == wb)
                            {
                                uint64_t mask = (bb < 63) ? ~((1ULL << (bb + 1)) - 1ULL) : 0ULL;
                                c &= mask;
                            }
                            while (c)
                            {
                                unsigned b = (unsigned)__builtin_ctzll(c);
                                int k = (w << 6) + (int)b;
                                if (k < N) { cf << i0 << '\t' << j0 << '\t' << k << '\n'; }
                                c &= (c - 1ULL);
                            }
                        }
                    }
                }
            }
            else if (cluster_size == 2)
            {
                cf << "# v0\tv1\n";
                for (int i0 = 0; i0 < N; ++i0)
                {
                    for (int j0 = i0 + 1; j0 < N; ++j0)
                    {
                        int wb = j0 >> 6; int bb = j0 & 63;
                        if (NB[(size_t)i0 * (size_t)W + (size_t)wb] & (1ULL << bb))
                            cf << i0 << '\t' << j0 << '\n';
                    }
                }
            }
            else
            {
                // k>=4 の列挙（一般 NB）
                cf << "# v0\tv1\t... (k=" << cluster_size << ")\n";
                dump_k_cliques_bitadj(cluster_size, NB, N, W, cf);
            }
            cf.close();
        }
        catch (...) { }
    }
}
// --- end worker refactor helpers ---

// ========================= simulate helpers (top-level) =========================
static uint64_t seed_from_obj(py::object seed)
{
    if (seed.is_none())
    {
        std::random_device rd;
        uint64_t t = (uint64_t)std::chrono::steady_clock::now().time_since_epoch().count();
        return ((uint64_t)rd() << 32) ^ (uint64_t)rd() ^ t;
    }
    else
    {
        return seed.cast<uint64_t>();
    }
}

static void parse_bins_option_or_default(py::object bins, int &bins_ra, int &bins_dec)
{
    bins_ra = 72;
    bins_dec = 36;
    if (bins.is_none())
    {
        return;
    }
    if (py::isinstance<py::int_>(bins))
    {
        int b = bins.cast<int>();
        if (b <= 0)
        {
            throw std::runtime_error("bins must be positive");
        }
        bins_ra = b;
        bins_dec = b;
        return;
    }
    if (py::isinstance<py::sequence>(bins))
    {
        py::sequence bs = bins;
        if (py::len(bs) != 2)
        {
            throw std::runtime_error("bins must be (ra, dec)");
        }
        bins_ra = py::cast<int>(bs[0]);
        bins_dec = py::cast<int>(bs[1]);
        if (bins_ra <= 0 || bins_dec <= 0)
        {
            throw std::runtime_error("bins must be positive");
        }
        return;
    }
    throw std::runtime_error("bins must be int or sequence of length 2");
}


static int decide_threads_param(py::object threads, bool debug, int M)
{
    int T = 1;
    if (!debug)
    {
        if (!threads.is_none())
        {
            int th = threads.cast<int>();
            if (th > 0)
            {
                T = th;
            }
        }
        else
        {
            unsigned hc = std::thread::hardware_concurrency();
            if (hc > 0)
            {
                T = (int)hc;
            }
        }
    }
    if (T < 1)
    {
        T = 1;
    }
    if (T > M)
    {
        T = M;
    }
    return T;
}


static void run_simulation_and_aggregate(
    int M,
    const std::vector<SiteC> &S,
    int N,
    double cos_thr,
    int cluster_size,
    double ra_lo,
    double ra_hi,
    int bins_ra,
    int bins_dec,
    bool return_histograms,
    py::array_t<int32_t> &counts,
    std::vector<long long> &H_evt2d,
    std::vector<long long> &H_tri2d,
    const std::function<void(double, double, int &, int &)> &bin2d,
    bool debug,
    bool progress,
    py::object threads,
    uint64_t seed_val,
    const std::string &log_text_dir)
{
    py::gil_scoped_release rel;

    const bool use_progress = (!debug) && progress && (M > 0);
    const int step = std::max(1, M / 200);
    auto start_steady = std::chrono::steady_clock::now();
    std::atomic<int> done(0);
    std::atomic<bool> cancel(false);

    std::vector<int32_t> counts_vec(M, 0);
    int T = decide_threads_param(threads, debug, M);

    std::mutex print_mtx;
    auto last_progress_time = start_steady;

    std::vector<std::vector<long long>> Evt2d_loc(T), Tri2d_loc(T);

    auto worker = [&](int tid, int t_begin, int t_end)
    {
        RNG rng(seed_val ^ (0x9E3779B97F4A7C15ULL * (uint64_t)(tid + 1)));
        ThreadBuf B(N);

        if (return_histograms)
        {
            Evt2d_loc[tid].assign((size_t)bins_ra * (size_t)bins_dec, 0);
            Tri2d_loc[tid].assign((size_t)bins_ra * (size_t)bins_dec, 0);
        }

        for (int t = t_begin; t < t_end; ++t)
        {
            if (cancel.load(std::memory_order_relaxed))
            {
                break;
            }
            if (debug)
            {
                std::fprintf(stderr, "t=%d/%d generate events (N=%d, sites=%zu)\n",
                             t + 1, M, N, S.size());
            }

            generate_events_and_cache(S, rng, B);

            update_event_histograms_if_needed(
                N, B.ra_deg_arr, B.dec_deg_arr, tid,
                return_histograms,
                bins_ra,
                bin2d, Evt2d_loc);

            if (N <= 128)
            {
                process_smallN_trial(t, M, N, tid, cluster_size, cos_thr,
                                     bins_ra, debug,
                                     return_histograms,
                                     bin2d, log_text_dir,
                                     B, Evt2d_loc,
                                     Tri2d_loc, counts_vec);
            }
            else
            {
                process_largeN_trial(t, M, N, tid, cluster_size, cos_thr,
                                     bins_ra, debug,
                                     return_histograms,
                                     bin2d, log_text_dir,
                                     B, Evt2d_loc,
                                     Tri2d_loc, counts_vec);
            }

            int c = done.fetch_add(1) + 1;
            update_progress_bar(use_progress, c, M, start_steady, last_progress_time, print_mtx, step);
        }
    };

    std::vector<std::thread> ths;
    ths.reserve(T);
    for (int tid = 0; tid < T; ++tid)
    {
        int t_begin = (int)((int64_t)M * tid / T);
        int t_end = (int)((int64_t)M * (tid + 1) / T);
        ths.emplace_back(worker, tid, t_begin, t_end);
    }

    while (!cancel.load(std::memory_order_relaxed) && done.load(std::memory_order_relaxed) < M)
    {
        {
            py::gil_scoped_acquire gil;
            if (PyErr_CheckSignals() != 0)
            {
                cancel.store(true, std::memory_order_relaxed);
                break;
            }
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(25));
    }

    for (auto &th : ths)
    {
        th.join();
    }

    if (cancel.load(std::memory_order_relaxed))
    {
        if (use_progress)
        {
            std::fputc('\n', stderr);
            std::fflush(stderr);
        }
        py::gil_scoped_acquire gil;
        PyErr_SetNone(PyExc_KeyboardInterrupt);
        throw py::error_already_set();
    }

    if (return_histograms)
    {
        size_t SZ = (size_t)bins_ra * (size_t)bins_dec;
        for (size_t tid = 0; tid < Evt2d_loc.size(); ++tid)
        {
            for (size_t i = 0; i < SZ; ++i) H_evt2d[i] += Evt2d_loc[tid][i];
            for (size_t i = 0; i < SZ; ++i) H_tri2d[i] += Tri2d_loc[tid][i];
        }
    }

    if (use_progress)
    {
        std::fputc('\n', stderr);
        std::fflush(stderr);
    }

    {
        py::gil_scoped_acquire gil;
        auto counts_v = counts.mutable_unchecked<1>();
        for (int t = 0; t < M; ++t)
        {
            counts_v(t) = counts_vec[t];
        }
    }
}
// ======================= end simulate helpers (top-level) =======================

py::dict simulate(
    int M,
    py::sequence sites, // [{lat, lon, zmax, n}, ...]
    double d_deg,
    py::object seed = py::none(), // 省略可：None→自動シード
    py::object bins = py::none(), // (ra, dec) または int または None
    bool return_histograms = false,
    py::object threads = py::none(),
    bool progress = true,
    bool debug = false,
    int cluster_size = 3)
{
    // 入力検証と共通前処理（しきい値や乱数シード、ビン設定など）を行い、
    // マルチスレッドで試行を回して集計した結果を Python dict で返します。
    // サイト配列を内部表現へ変換し、イベント総数 N を算出
    auto S = parse_sites(sites);
    const int N = total_events(S);
    if (N <= 0)
        throw std::runtime_error("Total events must be > 0");
    // 角距離しきい値の内積側基準（浮動誤差対策で微小緩和）
    const double cos_thr = std::cos(rad(d_deg)) - 1e-12;
    if (cluster_size < 2)
        throw std::runtime_error("cluster_size must be >= 2");

    uint64_t seed_val = seed_from_obj(seed);

    int bins_ra = 72, bins_dec = 36;
    parse_bins_option_or_default(bins, bins_ra, bins_dec);

    // 共通の bin 幅（RA は [-180,180) を採用。縫い目上の重複を避ける）
    const double ra_lo = -180.0;
    const double ra_hi = 180.0;
    const double wra = (ra_hi - ra_lo) / (double)bins_ra; // = 360/bins_ra
    const double wde = 180.0 / (double)bins_dec;

    std::string log_dir_path, log_text_dir, log_fig_dir;
    prepare_debug_dirs(debug, log_dir_path, log_text_dir, log_fig_dir);

    // 出力：トリプレット/クリーク数（counts 配列）。debug 時は M を 1000 にクランプ
    if (debug && M > 1000)
    {
        std::fprintf(stderr, "[debug] M clamped to 1000\n");
        M = 1000;
    }
    py::array_t<int32_t> counts(M);

    // 2D ヒスト（RA×Dec）。必要なときだけ確保
    std::vector<long long> H_evt2d, H_tri2d;
    if (return_histograms)
    {
        H_evt2d.assign((size_t)bins_ra * (size_t)bins_dec, 0);
        H_tri2d.assign((size_t)bins_ra * (size_t)bins_dec, 0);
    }

    auto bin2d = make_bin2d(bins_ra, bins_dec, ra_lo, ra_hi);

    run_simulation_and_aggregate(
        M, S, N, cos_thr, cluster_size,
        ra_lo, ra_hi, bins_ra, bins_dec,
        return_histograms,
        counts,
        H_evt2d, H_tri2d,
        bin2d,
        debug, progress, threads, seed_val,
        log_text_dir);

    // 返却（counts は常に返す。ヒストは要求時のみ）
    py::dict out;
    out["counts"] = counts;
    if (debug)
    {
        try
        {
            out["log_dir"] = py::str(log_dir_path);
        }
        catch (...)
        {
        }
    }

    // 共通エッジ（1D/2D で同じものを使う）
    if (return_histograms)
    {
        py::array_t<double> ra_edges(bins_ra + 1), dec_edges(bins_dec + 1);
        for (int i = 0; i <= bins_ra; ++i)
            ra_edges.mutable_at(i) = ra_lo + wra * (double)i; // [-180, 180]（右端は境界）
        for (int i = 0; i <= bins_dec; ++i)
            dec_edges.mutable_at(i) = -90.0 + wde * (double)i;
        out["ra_edges"] = ra_edges;
        out["dec_edges"] = dec_edges;
        out["ra_origin_deg"] = py::float_(ra_lo);
    }

    if (return_histograms)
    {
        py::array_t<long long> E({bins_dec, bins_ra});
        py::array_t<long long> T({bins_dec, bins_ra});
        auto Ev = E.mutable_unchecked<2>(), Tv = T.mutable_unchecked<2>();
        for (int id = 0; id < bins_dec; ++id)
            for (int ir = 0; ir < bins_ra; ++ir)
            {
                size_t idx = (size_t)id * (size_t)bins_ra + (size_t)ir;
                Ev(id, ir) = H_evt2d[idx];
                Tv(id, ir) = H_tri2d[idx];
            }
        out["events_hist2d"] = E;
        out["triplets_hist2d"] = T;
    }
    return out;
}

PYBIND11_MODULE(_core, m)
{
    m.doc() = "Triplet simulation: counts only by default; optional 1D/2D RA/Dec hist with unified bins; seed optional";
    m.def("simulate", &simulate,
          py::arg("M"),
          py::arg("sites"),
          py::arg("d_deg"),
          py::arg("seed") = py::none(),
          py::arg("bins") = py::none(),
          py::arg("return_histograms") = false,
          py::arg("threads") = py::none(),
          py::arg("progress") = true,
          py::arg("debug") = false,
          py::arg("cluster_size") = 3);
}
