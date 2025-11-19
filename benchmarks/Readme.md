- benchmarks/src/ → contains code to perform benchmarks
- benchmarks/results/ → cleaned, processed, human-readable results
- benchmarks/runs/ → raw execution output
- benchmarks/configs/ → reproducible configuration files


```
benchmarks/
│
├── src/                     # Benchmark source code (C++, Python, scripts)
│    ├── microbench/        # e.g. single-kernel timing
│    ├── macrobench/        # full-program benchmark drivers
│    └── utils/             # timing wrappers, measurement helpers
│
├── configs/                # benchmark-specific config templates
│
├── results/                # processed benchmark results (CSV, PNG)
│    ├── strong_scaling/
│    ├── weak_scaling/
│    ├── comm_profile/
│    └── sfc_load_balance/
│
└── runs/                   # raw outputs, per-job logs (gitignored)
└── <job_id>/
├── out.log
├── err.log
├── perf_raw.json
└── timings_raw.txt
```