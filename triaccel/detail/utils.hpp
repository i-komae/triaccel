// utils.hpp
//
// triaccel の C++ 実装で用いる共通ユーティリティを集約したヘッダ。
// 目的：
//  - ばらけやすい小粒な関数群（数学/ビット/ヒスト/ログ）を一箇所にまとめ、
//    _core.cpp から本筋でない処理を切り出して可読性を上げる。
//  - ファイル数を増やし過ぎないよう、細分化していたヘッダを統合。

#pragma once

#include <cstdint>
#include <cmath>
#include <string>
#include <cstdio>
#include <ctime>
#include <vector>
#include <mutex>
#include <chrono>
#include <functional>

// ============================= 数学ユーティリティ =============================
constexpr double PI = 3.14159265358979323846;
static inline double rad(double d) { return d * PI / 180.0; }
static inline double deg(double r) { return r * 180.0 / PI; }

// [-1,1] にクランプ
static inline double clamp_unit(double v)
{
    if (v < -1.0) return -1.0;
    if (v > 1.0) return 1.0;
    return v;
}

// 直交座標 → (RA, Dec) [度]
static inline void compute_ra_dec_deg_from_xyz(double x, double y, double z,
                                               double &ra_deg, double &dec_deg)
{
    double ra = std::atan2(y, x);
    if (ra < 0.0) ra += 2 * PI;
    double dec = std::asin(clamp_unit(z));
    ra_deg = deg(ra);
    dec_deg = deg(dec);
}

// RA を [-180,180) に折り返し
static inline double wrap_deg180(double ra_deg)
{
    double t = std::fmod(ra_deg + 180.0, 360.0);
    if (t < 0.0) t += 360.0;
    return t - 180.0;
}

// ヒスト用 RA 正規化（右端を含まないよう nextafter で微調整）
static inline double wrap_ra_hist(double ra_deg, double ra_hi)
{
    double ra_wrapped = wrap_deg180(ra_deg);
    if (!(ra_wrapped < ra_hi))
    {
        ra_wrapped = std::nextafter(ra_hi, 0.0);
    }
    return ra_wrapped;
}

// ============================= ビット演算ユーティリティ =============================
#if defined(_MSC_VER)
#include <intrin.h>
static inline int pop64(unsigned long long x) { return (int)__popcnt64(x); }
static inline int ctz64(unsigned long long x)
{
    unsigned long idx;
    _BitScanForward64(&idx, x);
    return (int)idx;
}
#else
static inline int pop64(uint64_t x) { return __builtin_popcountll(x); }
static inline int ctz64(uint64_t x) { return __builtin_ctzll(x); }
#endif

// ============================= デバッグ/ログ ユーティリティ =============================
#if defined(_WIN32)
#include <windows.h>
#include <direct.h>
#else
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/ioctl.h>
#include <unistd.h>
#endif

// 端末幅（stderr）を取得
static inline int get_terminal_width()
{
#if defined(_WIN32)
    CONSOLE_SCREEN_BUFFER_INFO csbi;
    if (GetConsoleScreenBufferInfo(GetStdHandle(STD_ERROR_HANDLE), &csbi))
        return (int)(csbi.srWindow.Right - csbi.srWindow.Left + 1);
    return 80;
#else
    struct winsize ws;
    if (ioctl(STDERR_FILENO, TIOCGWINSZ, &ws) == 0 && ws.ws_col > 0)
        return (int)ws.ws_col;
    return 80;
#endif
}

// 経過秒を HH:MM:SS 形式に整形
static inline std::string fmt_hms(double seconds)
{
    if (seconds < 0) seconds = 0;
    long long s = (long long)(seconds + 0.5);
    int h = (int)(s / 3600);
    int m = (int)((s % 3600) / 60);
    int sec = (int)(s % 60);
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%02d:%02d:%02d", h, m, sec);
    return std::string(buf);
}

// 簡易ディレクトリ作成（'/' 区切り前提）
static inline void ensure_dirs(const std::string &path)
{
    if (path.empty()) return;
    std::string tmp;
    tmp.reserve(path.size());
    for (size_t i = 0; i < path.size(); ++i)
    {
        char c = path[i];
        tmp.push_back(c);
        if (c == '/' || i + 1 == path.size())
        {
            if (tmp.size() == 1 && tmp[0] == '/') continue;
#if defined(_WIN32)
            int rc = _mkdir(tmp.c_str());
#else
            int rc = mkdir(tmp.c_str(), 0755);
#endif
            (void)rc;
        }
    }
}

// UTC の ISO8601 風タイムスタンプ（ファイル名用に ':' → '-' など）
static inline std::string iso8601_utc_now()
{
    std::time_t t = std::time(nullptr);
    std::tm tm{};
#if defined(_WIN32)
    gmtime_s(&tm, &t);
#else
    gmtime_r(&t, &tm);
#endif
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%04d-%02d-%02dT%02d-%02d-%02dZ",
                  tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday,
                  tm.tm_hour, tm.tm_min, tm.tm_sec);
    return std::string(buf);
}

// 進捗バーを標準エラーに表示（高頻度更新を抑制）
static inline void update_progress_bar(
    bool use_progress,
    int c,
    int M,
    const std::chrono::steady_clock::time_point &start_steady,
    std::chrono::steady_clock::time_point &last_progress_time,
    std::mutex &print_mtx,
    int step)
{
    if (!use_progress) return;
    auto now_steady = std::chrono::steady_clock::now();
    std::lock_guard<std::mutex> lk(print_mtx);
    auto since = std::chrono::duration_cast<std::chrono::seconds>(now_steady - last_progress_time).count();
    if ((c % step) != 0 && c != M && since < 1) return;
    last_progress_time = now_steady;

    double pct = (100.0 * c) / (double)M;
    double elapsed_s = std::chrono::duration_cast<std::chrono::seconds>(now_steady - start_steady).count();
    double rem_s = (c > 0) ? (elapsed_s * ((double)M / (double)c - 1.0)) : 0.0;
    if (rem_s < 0) rem_s = 0.0;

    std::string meta = [&]{
        std::string elapsed_str = fmt_hms(elapsed_s);
        std::string eta_str = fmt_hms(rem_s);
        char meta_buf[128];
        std::snprintf(meta_buf, sizeof(meta_buf), "%6.2f%% (%d/%d)  elapsed %s  ETA %s",
                      pct, c, M, elapsed_str.c_str(), eta_str.c_str());
        return std::string(meta_buf);
    }();

    int cols = get_terminal_width();
    int width = cols - 4 - (int)meta.size();
    if (width < 10) width = 10;
    if (width > 200) width = 200;

    int filled = (int)((pct / 100.0) * width);
    if (filled < 0) filled = 0;
    if (filled > width) filled = width;

    std::fprintf(stderr, "\r[");
    for (int ii = 0; ii < filled; ++ii) std::fputc('#', stderr);
    for (int ii = filled; ii < width; ++ii) std::fputc('.', stderr);
    std::fprintf(stderr, "] %s ", meta.c_str());
    std::fflush(stderr);
}

// デバッグ用ログディレクトリの準備（例: log/2024-.../{text,fig}）
static inline void prepare_debug_dirs(bool debug,
                                      std::string &log_dir_path,
                                      std::string &log_text_dir,
                                      std::string &log_fig_dir)
{
    log_dir_path.clear();
    log_text_dir.clear();
    log_fig_dir.clear();
    if (!debug) return;
    try
    {
        auto ts = iso8601_utc_now();
        log_dir_path = std::string("log/") + ts;
        log_text_dir = log_dir_path + "/text";
        log_fig_dir = log_dir_path + "/fig";
        ensure_dirs(log_text_dir);
        ensure_dirs(log_fig_dir);
    }
    catch (...) { }
}

// ============================= ヒストグラムユーティリティ =============================
// 角度値を [lo, hi) 等間隔ビンにカウントアップ
static inline void bump_hist_1d(double val_deg, double lo, double hi, int bins, std::vector<long long> &H)
{
    double w = (hi - lo) / bins;
    int i = (int)std::floor((val_deg - lo) / w);
    if (i < 0) i = 0;
    if (i >= bins) i = bins - 1;
    H[(size_t)i] += 1;
}

// RA/Dec の 2D ビン関数（RA は [-180, 180)）
static inline std::function<void(double, double, int &, int &)> make_bin2d(int bins_ra,
                                                                           int bins_dec,
                                                                           double ra_lo,
                                                                           double ra_hi)
{
    double wra = (ra_hi - ra_lo) / (double)bins_ra;
    double wde = 180.0 / (double)bins_dec;
    return [=](double ra_deg, double dec_deg, int &ir, int &id)
    {
        double ra_wrapped = wrap_ra_hist(ra_deg, ra_hi);
        ir = (int)std::floor((ra_wrapped - ra_lo) / wra);
        if (ir < 0) ir = 0;
        if (ir >= bins_ra) ir = bins_ra - 1;
        int id_tmp = (int)std::floor((dec_deg + 90.0) / wde);
        if (id_tmp < 0) id_tmp = 0;
        if (id_tmp >= bins_dec) id_tmp = bins_dec - 1;
        id = id_tmp;
    };
}
