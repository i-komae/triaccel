
import time
import triaccel
import numpy as np

def benchmark():
    # Setup parameters similar to README example but larger M for measurable time
    sites = [
        {"lat": 39.3, "lon": -112.9, "zmax": 55.0, "n": 64},
        {"lat": -35.2, "lon": -69.2, "zmax": 60.0, "n": 64},
    ]
    # Total N = 128 (fits in smallN optimization path)

    M = 1_000_000
    d_deg = 5.0

    print(f"Benchmarking with M={M}, N=128 (smallN path)...")
    start = time.time()
    triaccel.simulate(
        M=M,
        sites=sites,
        d_deg=d_deg,
        bins=None, # No histogram for pure core speed test
        return_histograms=False,
        cluster_size=3,
        progress=False
    )
    end = time.time()
    print(f"Time (N=128): {end - start:.4f} seconds")

    # Larger N case
    sites_large = [
        {"lat": 39.3, "lon": -112.9, "zmax": 55.0, "n": 100},
        {"lat": -35.2, "lon": -69.2, "zmax": 60.0, "n": 100},
    ]
    # Total N = 200 (largeN path)

    print(f"Benchmarking with M={M}, N=200 (largeN path)...")
    start = time.time()
    triaccel.simulate(
        M=M,
        sites=sites_large,
        d_deg=d_deg,
        bins=None,
        return_histograms=False,
        cluster_size=3,
        progress=False
    )
    end = time.time()
    print(f"Time (N=200): {end - start:.4f} seconds")

if __name__ == "__main__":
    benchmark()
