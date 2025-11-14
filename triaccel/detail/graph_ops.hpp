// graph_ops.hpp
//
// 隣接行列の構築（小規模=2語 / 大規模=任意幅）と、
// ペア・三角形・一般 k-クリークのカウント／列挙ヘルパ群を集約。
// 役割を _core.cpp から切り出して見通しを良くする。

#pragma once

#include <vector>
#include <cstdint>
#include <functional>
#include <cmath>
#include "utils.hpp" // pop64, ctz64, clamp, wrap_ra_hist, bump_hist_1d など

// --------------------------- 小規模 N (<=128) ---------------------------
inline void build_masks_u128(const std::vector<double> &x,
                             const std::vector<double> &y,
                             const std::vector<double> &z,
                             double cos_thr,
                             std::vector<uint64_t> &m0,
                             std::vector<uint64_t> &m1)
{
    int N = (int)x.size();
    for (int i = 0; i < N; ++i) { m0[i] = m1[i] = 0ULL; }
    for (int i = 0; i < N - 1; ++i)
    {
        double xi = x[i], yi = y[i], zi = z[i];
        for (int j = i + 1; j < N; ++j)
        {
            double dot = xi * x[j] + yi * y[j] + zi * z[j];
            if (dot >= cos_thr)
            {
                if (j < 64) m0[i] |= (1ULL << j); else m1[i] |= (1ULL << (j - 64));
                if (i < 64) m0[j] |= (1ULL << i); else m1[j] |= (1ULL << (i - 64));
            }
        }
    }
}

inline long long count_triangles_u128(const std::vector<uint64_t> &m0,
                                      const std::vector<uint64_t> &m1)
{
    int N = (int)m0.size();
    long long cnt = 0;
    for (int i = 0; i < N - 2; ++i)
    {
        for (int j = i + 1; j < N - 1; ++j)
        {
            bool ij = (j < 64) ? (((m0[i] >> j) & 1ULL) != 0ULL)
                               : (((m1[i] >> (j - 64)) & 1ULL) != 0ULL);
            if (!ij) continue;
            uint64_t c0 = m0[i] & m0[j];
            uint64_t c1 = m1[i] & m1[j];
            if (j < 63)
                c0 &= ~(((1ULL) << (j + 1)) - 1ULL);
            else if (j == 63)
                c0 = 0ULL;
            else
            {
                c0 = 0ULL; int jb = j - 64; c1 &= ~(((1ULL) << (jb + 1)) - 1ULL);
            }
            cnt += pop64(c0) + pop64(c1);
        }
    }
    return cnt;
}

inline long long count_pairs_u128(const std::vector<uint64_t> &m0,
                                  const std::vector<uint64_t> &m1,
                                  int N,
                                  std::vector<char> *in_cluster,
                                  bool debug)
{
    long long deg_sum = 0;
    for (int i = 0; i < N; ++i)
    {
        int di = pop64(m0[i]) + pop64(m1[i]);
        deg_sum += di;
        if (debug && in_cluster && di > 0) (*in_cluster)[i] = 1;
    }
    return deg_sum / 2;
}

inline void enumerate_triangles_hist_from_masks_u128(
    int N, int tid,
    const std::vector<double> &x,
    const std::vector<double> &y,
    const std::vector<double> &z,
    const std::vector<uint64_t> &m0,
    const std::vector<uint64_t> &m1,
    bool mark_members,
    std::vector<char> &in_cluster,
    bool return_hist2d,
    int bins_ra,
    const std::function<void(double, double, int &, int &)> &bin2d,
    std::vector<std::vector<long long>> &Tri2d_loc,
    bool tri2d_push_centroid)
{
    for (int i = 0; i < N - 2; ++i)
    {
        for (int j = i + 1; j < N - 1; ++j)
        {
            bool ij = (j < 64) ? (((m0[i] >> j) & 1ULL) != 0ULL)
                               : (((m1[i] >> (j - 64)) & 1ULL) != 0ULL);
            if (!ij) continue;
            uint64_t c0 = m0[i] & m0[j];
            uint64_t c1 = m1[i] & m1[j];
            if (j < 63) c0 &= ~(((1ULL) << (j + 1)) - 1ULL);
            else if (j == 63) c0 = 0ULL;
            else { c0 = 0ULL; int jb = j - 64; c1 &= ~(((1ULL) << (jb + 1)) - 1ULL); }
            while (c0)
            {
                int k = ctz64(c0); c0 &= (c0 - 1);
                if (mark_members) { in_cluster[i] = in_cluster[j] = in_cluster[k] = 1; }
                // push histogram
                double vx = x[i] + x[j] + x[k];
                double vy = y[i] + y[j] + y[k];
                double vz = z[i] + z[j] + z[k];
                double r = std::sqrt(vx*vx + vy*vy + vz*vz); if (r == 0) r = 1.0;
                vx/=r; vy/=r; vz/=r;
                double ra_deg, dec_deg; compute_ra_dec_deg_from_xyz(vx, vy, vz, ra_deg, dec_deg);
                if (return_hist2d && tri2d_push_centroid)
                {
                    int ir, id; bin2d(ra_deg, dec_deg, ir, id);
                    Tri2d_loc[tid][(size_t)id * (size_t)bins_ra + (size_t)ir] += 1;
                }
            }
            const int base = 64;
            while (c1)
            {
                int kb = ctz64(c1); c1 &= (c1 - 1); int k = base + kb; if (k >= N) continue;
                if (mark_members) { in_cluster[i] = in_cluster[j] = in_cluster[k] = 1; }
                double vx = x[i] + x[j] + x[k];
                double vy = y[i] + y[j] + y[k];
                double vz = z[i] + z[j] + z[k];
                double r = std::sqrt(vx*vx + vy*vy + vz*vz); if (r == 0) r = 1.0;
                vx/=r; vy/=r; vz/=r;
                double ra_deg, dec_deg; compute_ra_dec_deg_from_xyz(vx, vy, vz, ra_deg, dec_deg);
                if (return_hist2d && tri2d_push_centroid)
                {
                    int ir, id; bin2d(ra_deg, dec_deg, ir, id);
                    Tri2d_loc[tid][(size_t)id * (size_t)bins_ra + (size_t)ir] += 1;
                }
            }
        }
    }
}

inline long long count_k_cliques_u128(
    int K,
    const std::vector<uint64_t> &m0,
    const std::vector<uint64_t> &m1,
    int N,
    std::vector<char> *in_cluster,
    bool debug)
{
    long long kcnt = 0;
    std::vector<int> R; R.reserve(K);
    auto mask_gt = [&](uint64_t &a0, uint64_t &a1, int thr)
    {
        if (thr < 63) a0 &= ~(((1ULL) << (thr + 1)) - 1ULL);
        else if (thr == 63) { a0 = 0ULL; }
        else { a0 = 0ULL; int jb = thr - 64; if (jb < 63) a1 &= ~(((1ULL) << (jb + 1)) - 1ULL); else a1 = 0ULL; }
    };
    std::function<void(uint64_t, uint64_t)> dfs = [&](uint64_t P0, uint64_t P1)
    {
        int rsz = (int)R.size(); int psize = pop64(P0) + pop64(P1);
        if (rsz + psize < K) return;
        if (rsz == K)
        {
            kcnt += 1; if (debug && in_cluster) for (int v : R) (*in_cluster)[v] = 1; return;
        }
        uint64_t Q0 = P0, Q1 = P1;
        while (Q0 || Q1)
        {
            int v; if (Q0) { int b = ctz64(Q0); Q0 &= (Q0 - 1); v = b; }
            else { int b = ctz64(Q1); Q1 &= (Q1 - 1); v = 64 + b; if (v >= N) continue; }
            R.push_back(v);
            uint64_t NP0 = m0[v] & P0; uint64_t NP1 = m1[v] & P1; mask_gt(NP0, NP1, v);
            dfs(NP0, NP1); R.pop_back();
        }
    };
    for (int v0 = 0; v0 < N; ++v0)
    {
        uint64_t P0 = m0[v0], P1 = m1[v0]; mask_gt(P0, P1, v0);
        R.push_back(v0); dfs(P0, P1); R.pop_back();
    }
    return kcnt;
}

inline void enumerate_pairs_and_push_hist_u128(
    int N, int tid,
    const std::vector<uint64_t> &m0,
    const std::vector<uint64_t> &m1,
    const std::vector<double> &x,
    const std::vector<double> &y,
    const std::vector<double> &z,
    bool return_hist2d,
    int bins_ra,
    const std::function<void(double, double, int &, int &)> &bin2d,
    std::vector<std::vector<long long>> &H_2d_loc)
{
    for (int i = 0; i < N - 1; ++i)
    {
        uint64_t a0 = (i < 63) ? (m0[i] & ~(((1ULL) << (i + 1)) - 1ULL)) : 0ULL;
        while (a0)
        {
            int j = ctz64(a0); a0 &= (a0 - 1);
            std::vector<int> R = {i, j};
            double vx = x[i] + x[j]; double vy = y[i] + y[j]; double vz = z[i] + z[j];
            double r = std::sqrt(vx*vx + vy*vy + vz*vz); if (r == 0) r = 1.0; vx/=r; vy/=r; vz/=r;
            double ra_deg, dec_deg; compute_ra_dec_deg_from_xyz(vx, vy, vz, ra_deg, dec_deg);
            if (return_hist2d) { int ir, id; bin2d(ra_deg, dec_deg, ir, id);
                H_2d_loc[tid][(size_t)id * (size_t)bins_ra + (size_t)ir] += 1; }
        }
        uint64_t a1 = m1[i]; if (i >= 64) { int ib = i - 64; if (ib < 63) a1 &= ~(((1ULL) << (ib + 1)) - 1ULL); else a1 = 0ULL; }
        const int base = 64;
        while (a1)
        {
            int jb = ctz64(a1); a1 &= (a1 - 1); int j = base + jb; if (j >= N) continue;
            double vx = x[i] + x[j]; double vy = y[i] + y[j]; double vz = z[i] + z[j];
            double r = std::sqrt(vx*vx + vy*vy + vz*vz); if (r == 0) r = 1.0; vx/=r; vy/=r; vz/=r;
            double ra_deg, dec_deg; compute_ra_dec_deg_from_xyz(vx, vy, vz, ra_deg, dec_deg);
            if (return_hist2d) { int ir, id; bin2d(ra_deg, dec_deg, ir, id);
                H_2d_loc[tid][(size_t)id * (size_t)bins_ra + (size_t)ir] += 1; }
        }
    }
}

// --------------------------- 任意 N（ビット隣接） ---------------------------
inline void masks_to_NB(const std::vector<uint64_t> &m0,
                        const std::vector<uint64_t> &m1,
                        int N,
                        std::vector<uint64_t> &NB,
                        int &W)
{ W = (N + 63) / 64; NB.assign((size_t)N * (size_t)W, 0ULL);
  for (int i = 0; i < N; ++i) { if (W >= 1) NB[(size_t)i * (size_t)W + 0] = m0[i]; if (W >= 2) NB[(size_t)i * (size_t)W + 1] = m1[i]; } }

inline bool bit_test_edge(const std::vector<uint64_t> &NB, int N, int W, int a, int b)
{ (void)N; int wb = b >> 6; int bb = b & 63; return ((NB[(size_t)a * (size_t)W + (size_t)wb] >> bb) & 1ULL) != 0ULL; }

inline long long count_pairs_bitadj(const std::vector<uint64_t> &NB, int N, int W,
                                    std::vector<char> *in_cluster, bool debug)
{
    long long pair_cnt = 0; for (int i = 0; i < N; ++i) for (int j = i + 1; j < N; ++j)
        if (bit_test_edge(NB, N, W, i, j)) { pair_cnt += 1; if (debug && in_cluster) { (*in_cluster)[i] = 1; (*in_cluster)[j] = 1; } }
    return pair_cnt;
}

inline long long count_and_hist_triangles_bitadj(
    const std::vector<uint64_t> &NB, int N, int W, int tid,
    const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &z,
    bool debug, std::vector<char> &in_cluster,
    bool return_hist2d,
    int bins_ra,
    const std::function<void(double, double, int &, int &)> &bin2d,
    std::vector<std::vector<long long>> &Tri2d_loc,
    bool tri2d_push_centroid)
{
    long long tri_cnt = 0;
    for (int i2 = 0; i2 < N - 2; ++i2)
    {
        for (int j2 = i2 + 1; j2 < N - 1; ++j2)
        {
            if (!bit_test_edge(NB, N, W, i2, j2)) continue;
            int jword = j2 >> 6; int jbit = j2 & 63;
            for (int w = jword; w < W; ++w)
            {
                uint64_t mask = NB[(size_t)i2 * (size_t)W + (size_t)w] & NB[(size_t)j2 * (size_t)W + (size_t)w];
                if (w == jword)
                {
                    uint64_t lowmask = (jbit == 63) ? ~0ULL : ((1ULL << (jbit + 1)) - 1ULL);
                    mask &= ~lowmask;
                }
                while (mask)
                {
                    int b2 = ctz64(mask); mask &= (mask - 1); int k2 = (w << 6) + b2; if (k2 >= N) continue;
                    tri_cnt += 1; if (debug) { in_cluster[i2] = in_cluster[j2] = in_cluster[k2] = 1; }
                    if (return_hist2d)
                    {
                        double vx = x[i2] + x[j2] + x[k2]; double vy = y[i2] + y[j2] + y[k2]; double vz = z[i2] + z[j2] + z[k2];
                        double r = std::sqrt(vx*vx + vy*vy + vz*vz); if (r == 0) r = 1.0; vx/=r; vy/=r; vz/=r;
                        double ra_deg, dec_deg; compute_ra_dec_deg_from_xyz(vx, vy, vz, ra_deg, dec_deg);
                        if (tri2d_push_centroid) { int ir, id; bin2d(ra_deg, dec_deg, ir, id);
                            Tri2d_loc[tid][(size_t)id * (size_t)bins_ra + (size_t)ir] += 1; }
                    }
                }
            }
        }
    }
    return tri_cnt;
}

inline long long count_k_cliques_bitadj(
    int K,
    const std::vector<uint64_t> &NB,
    int N,
    int W,
    std::vector<char> *in_cluster,
    bool debug,
    int tid,
    const std::vector<double> *x,
    const std::vector<double> *y,
    const std::vector<double> *z,
    bool return_hist2d,
    int bins_ra,
    const std::function<void(double, double, int &, int &)> *bin2d,
    std::vector<std::vector<long long>> *H_2d_loc)
{
    long long kcnt = 0; std::vector<int> R; R.reserve(K);
    auto pop_bits = [&](const std::vector<uint64_t> &A)->long long { long long s = 0; for (uint64_t w : A) s += pop64(w); return s; };
    auto mask_gt_assign = [&](std::vector<uint64_t> &P, int thr)
    {
        int wthr = thr >> 6; int bthr = thr & 63;
        for (int w = 0; w < wthr; ++w) P[(size_t)w] = 0ULL;
        uint64_t lowmask = (bthr == 63) ? ~0ULL : ((1ULL << (bthr + 1)) - 1ULL);
        P[(size_t)wthr] &= ~lowmask;
    };
    std::function<void(std::vector<uint64_t>&)> dfs;
    dfs = [&](std::vector<uint64_t> &P)
    {
        int rsz = (int)R.size(); long long psize = pop_bits(P);
        if (rsz + psize < K) return;
        if (rsz == K)
        {
            kcnt += 1; if (debug && in_cluster) for (int v : R) (*in_cluster)[v] = 1;
            if (return_hist2d)
            {
                double vx = 0.0, vy = 0.0, vz = 0.0; for (int v : R) { vx += (*x)[v]; vy += (*y)[v]; vz += (*z)[v]; }
                double r = std::sqrt(vx*vx + vy*vy + vz*vz); if (r == 0) r = 1.0; vx/=r; vy/=r; vz/=r;
                double ra_deg, dec_deg; compute_ra_dec_deg_from_xyz(vx, vy, vz, ra_deg, dec_deg);
                int ir, id; (*bin2d)(ra_deg, dec_deg, ir, id);
                (*H_2d_loc)[tid][(size_t)id * (size_t)bins_ra + (size_t)ir] += 1;
            }
            return;
        }
        std::vector<uint64_t> Q = P;
        for (int w = 0; w < W; ++w)
        {
            uint64_t mask = Q[(size_t)w];
            while (mask)
            {
                int b = ctz64(mask); mask &= (mask - 1); int v = (w << 6) + b; if (v >= N) continue;
                R.push_back(v);
                std::vector<uint64_t> NP(W);
                for (int ww = 0; ww < W; ++ww) NP[(size_t)ww] = P[(size_t)ww] & NB[(size_t)v * (size_t)W + (size_t)ww];
                mask_gt_assign(NP, v);
                dfs(NP);
                R.pop_back();
            }
        }
    };
    for (int v0 = 0; v0 < N; ++v0)
    {
        std::vector<uint64_t> P(W); for (int w = 0; w < W; ++w) P[(size_t)w] = NB[(size_t)v0 * (size_t)W + (size_t)w];
        mask_gt_assign(P, v0); R.push_back(v0); dfs(P); R.pop_back();
    }
    return kcnt;
}

inline void enumerate_pairs_and_push_hist_bitadj(
    const std::vector<uint64_t> &NB, int N, int W, int tid,
    const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &z,
    bool return_hist2d,
    int bins_ra,
    const std::function<void(double, double, int &, int &)> &bin2d,
    std::vector<std::vector<long long>> &H_2d_loc)
{
    for (int i = 0; i < N - 1; ++i)
        for (int j = i + 1; j < N; ++j)
            if (bit_test_edge(NB, N, W, i, j))
            {
                double vx = x[i] + x[j]; double vy = y[i] + y[j]; double vz = z[i] + z[j];
                double r = std::sqrt(vx*vx + vy*vy + vz*vz); if (r == 0) r = 1.0; vx/=r; vy/=r; vz/=r;
                double ra_deg, dec_deg; compute_ra_dec_deg_from_xyz(vx, vy, vz, ra_deg, dec_deg);
                if (return_hist2d) { int ir, id; bin2d(ra_deg, dec_deg, ir, id);
                    H_2d_loc[tid][(size_t)id * (size_t)bins_ra + (size_t)ir] += 1; }
            }
}

// end graph ops
