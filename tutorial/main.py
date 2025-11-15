import triaccel

N_TRIALS = 100_000_000
SAVE_HIST = False

sites = [{"lat": 39.3, "lon": -112.9, "zmax": 55.0, "n": 32},
         {"lat": -35.2, "lon": -69.2, "zmax": 60.0, "n": 35}]

res = triaccel.simulate(
    M=N_TRIALS,
    sites=sites,
    d_deg=5.0,
    bins=(144, 72),
    return_histograms=SAVE_HIST,
    debug=False,
    cluster_size=3,
)

# トリプレット確率分布を保存
triaccel.viz.write_triplet_prob_sigma(res, n=30, outdir="results", filename=f"counts_{N_TRIALS:.0e}.npy")

# 図を描く
# トリプレット数分布（log軸で保存）
triaccel.viz.count_hist(res, log=True, outdir="figs",
                        filename=f"triplet_counts_{N_TRIALS:.0e}")

if SAVE_HIST:
    # 1次元ヒスト（イベント/トリプレット）
    triaccel.viz.hist1d(res, what="events",   outdir="figs",
                        filename=f"events_1d_{N_TRIALS:.0e}")
    triaccel.viz.hist1d(res, what="triplets", outdir="figs",
                        filename=f"triplets_1d_{N_TRIALS:.0e}")
    # Hammer 図
    triaccel.viz.hammer(res, what="events",   outdir="figs",
                        filename=f"events_map_{N_TRIALS:.0e}")
    triaccel.viz.hammer(res, what="triplets", outdir="figs",
                        filename=f"triplets_map_{N_TRIALS:.0e}")
