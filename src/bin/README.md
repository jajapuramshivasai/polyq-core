# to profile code run
cargo flamegraph --bench sim_bench
### view target/flamegraph.svg
# to run benchmarks run
cargo bench
### view target/criterion/report/index.html

# current state

## run 3
 - rust QASM simulation time: 6.91563375s
 
 {
  "simulation_time_seconds": 7.423377990722656,
  "system_info": {
    "platform": "Darwin 25.2.0",
    "architecture": "arm64",
    "processor": "arm",
    "cpu_count": 8,
    "total_ram_gb": 8.0,
    "available_ram_gb": 1.31,
    "cpu_brand": "Apple M1",
    "cpu_hz": "N/A"
  },
  "date_time": "2026-03-08T00:02:13.099149"
}


