# to profile code run
cargo flamegraph --bench sim_bench
### view target/flamegraph.svg
# to run benchmarks run
cargo bench
### view target/criterion/report/index.html

# current state

## run 2
 - rust QASM simulation time: 121.274100167s
 
 {
  "simulation_time_seconds": 5.966247797012329,
  "system_info": {
    "platform": "Darwin 25.2.0",
    "architecture": "arm64",
    "processor": "arm",
    "cpu_count": 8,
    "total_ram_gb": 8.0,
    "available_ram_gb": 1.25,
    "cpu_brand": "Apple M1",
    "cpu_hz": "N/A"
  },
  "date_time": "2026-03-07T22:55:18.059837"
}


