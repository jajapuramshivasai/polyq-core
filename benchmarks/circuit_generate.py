from qiskit import QuantumCircuit
from qiskit_aer import AerSimulator
from qiskit.qasm2 import dumps
import numpy as np
import os
import time
import datetime
from collections import Counter
from qiskit.transpiler.passes import ResourceEstimation
from IPython.display import display
import platform
import psutil

def get_system_info():
    info = {
        "platform": f"{platform.system()} {platform.release()}",
        "architecture": platform.machine(),
        "processor": platform.processor(),
        "cpu_count": os.cpu_count(),
        "total_ram_gb": round(psutil.virtual_memory().total / (1024**3), 2),
        "available_ram_gb": round(psutil.virtual_memory().available / (1024**3), 2)
    }
    try:
        import cpuinfo
        cpu = cpuinfo.get_cpu_info()
        info["cpu_brand"] = cpu.get('brand_raw', 'N/A')
        info["cpu_hz"] = cpu.get('hz_advertised_friendly', 'N/A')
    except ImportError:
        info["cpu_brand"] = 'N/A'
        info["cpu_hz"] = 'N/A'
    return info

def generate_random_clifford_t_circuit(num_qubits, h_count, t_count, s_count=0, z_count=0, cz_count=0):
    """Generate a random quantum circuit using {H, Z, S, T, CZ} gates."""
    qc = QuantumCircuit(num_qubits)
    qc.h(0)
    for i in range(num_qubits-1):
        qc.cz(i, i + 1)


    # Add random H, S, Z gates
    for _ in range(h_count):
        qc.h(np.random.randint(0, num_qubits))
    for _ in range(s_count):
        qc.s(np.random.randint(0, num_qubits))
    for _ in range(z_count):
        qc.z(np.random.randint(0, num_qubits))

    # Add exactly t_count T gates, distributed randomly
    t_qubits = np.random.choice(num_qubits, size=t_count, replace=True)
    for qubit in t_qubits:
        qc.t(qubit)

    # Add random CZ gates
    for _ in range(cz_count):
        control = np.random.randint(0, num_qubits)
        target = np.random.randint(0, num_qubits)
        while target == control:
            target = np.random.randint(0, num_qubits)
        qc.cz(control, target)

    return qc


import shutil
import os

def clean_dataset_folders():
    # Remove all contents of 'dataset' and '../dataset' folders
    for folder in ['dataset', '../dataset']:
        if os.path.exists(folder):
            for entry in os.listdir(folder):
                entry_path = os.path.join(folder, entry)
                if os.path.isdir(entry_path):
                    shutil.rmtree(entry_path)
                else:
                    os.remove(entry_path)
    print("All contents of 'dataset' and '../dataset' have been deleted.")

def run_sample_1():
    # Set fixed RNG seed for reproducibility
    np.random.seed(42)

    # Parameters for random circuit
    num_qubits = 29
    h_count = 10
    t_count = 5
    s_count = 10
    z_count = 25
    cz_count = 10

    random_qc = generate_random_clifford_t_circuit(num_qubits, h_count, t_count, s_count, z_count, cz_count)

    # Save to 'test' directory with 'test' name
    test_dir = "test"
    os.makedirs(test_dir, exist_ok=True)
    filename_base = "test"
    d = test_dir
    os.makedirs(d, exist_ok=True)

    # Export to QASM2
    qasm_output = dumps(random_qc)
    qasm_filename = f"{d}/{filename_base}.qasm"
    with open(qasm_filename, 'w') as f:
        f.write(qasm_output)
    print(f"QASM saved to: {qasm_filename}")

    # Simulate and get statevector, benchmark simulation time
    simulator = AerSimulator(method='statevector')
    random_qc.save_statevector()
    start_sim = time.time()
    result = simulator.run(random_qc).result()
    elapsed_sim = time.time() - start_sim
    statevector = result.get_statevector()
    print(f"\nQiskit simulation time: {elapsed_sim:.4f} seconds")
    print(f"First element of statevector: {statevector[0]}")

    # Save simulation time and system info as JSON
    import json
    siminfo = {
        "simulation_time_seconds": elapsed_sim,
        "system_info": get_system_info(),
        "date_time": datetime.datetime.now().isoformat()
    }
    siminfo_filename = f"{d}/{filename_base}_siminfo.json"
    with open(siminfo_filename, 'w') as f:
        json.dump(siminfo, f, indent=2)
    print(f"Simulation info saved as: {siminfo_filename}")

    # Print gate counts
    gate_counts = Counter(random_qc.count_ops())
    print("Gate counts:")
    for gate in ['h', 't', 's', 'z', 'cz']:
        print(f"{gate.upper()}: {gate_counts.get(gate, 0)}")
    print(random_qc)
    
    
def generate_sweep_benchmark():
    # Set fixed RNG seed for reproducibility
    np.random.seed(42)

    # Define parameter ranges
    qubit_counts = [4, 8, 12, 16, 20, 24, 29]
    t_counts = [0, 5, 9, 11]

    # Create dataset directory
    dataset_dir = "dataset"
    os.makedirs(dataset_dir, exist_ok=True)

    # Loop over parameter combinations
    for num_qubits in qubit_counts:
        for t_count in t_counts:
            h_count = num_qubits // 3
            s_count = num_qubits // 3
            z_count = num_qubits // 3
            cz_count = num_qubits // 4

            random_qc = generate_random_clifford_t_circuit(num_qubits, h_count, t_count, s_count, z_count, cz_count)

            # Save to dataset directory with descriptive name
            filename_base = f"qc_{num_qubits}q_{t_count}T"
            d = dataset_dir
            os.makedirs(d, exist_ok=True)

            # Export to QASM2
            qasm_output = dumps(random_qc)
            qasm_filename = f"{d}/{filename_base}.qasm"
            with open(qasm_filename, 'w') as f:
                f.write(qasm_output)
            print(f"QASM saved to: {qasm_filename}")