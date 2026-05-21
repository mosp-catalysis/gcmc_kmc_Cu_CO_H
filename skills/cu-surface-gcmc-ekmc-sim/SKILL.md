---
name: cu-surface-gcmc-ekmc-sim
description: "Automates Cu surface GCMC+EKMC simulations. Clones repo, checks convergence, and dynamically limits runs by physical time. Robust pure-python post-processing."
metadata:
  openclaw:
    emoji: "🧪"
    requires:
      bins:
        - git
        - gfortran
        - make
        - awk
        - sed
        - python
        - nohup
---

# Cu Surface GCMC+KMC Full Workflow

This skill deploys a GCMC+EKMC simulation end-to-end. It enforces safety limits, handles dynamic surfaces, and executes a robust, pure-bash and pure-python workflow.

## 0. Target Parameters (Crucial Step)
Before generating any bash scripts, you **MUST** determine two variables from the user's prompt:
1. **`{{TARGET_SURFACE}}`**: Map the user's request to one of these strict canonical names: `Cu100`, `Cu111`, `Cu1111`, `Cu311`, `Cu511`, `Cu711`, `Cu911`. If it cannot be mapped, STOP and ask the user.
2. **`{{TARGET_TIME}}`**: The simulation target time in seconds. If the user doesn't provide it, STOP and ask. **NO DEFAULT ALLOWED.** Extract purely the numerical value.

*In all code blocks below, you MUST replace the literal strings `{{TARGET_SURFACE}}` and `{{TARGET_TIME}}` with the user's actual choices.*

## 1. Clone Repository & Compile
Clone from the strictly fixed URL and compile. Use the `shell` tool:

```bash
# MUST be a pure URL, do not use markdown link formatting
BASE_DIR="${HOME}/zsq/openclaw"
mkdir -p "$BASE_DIR"
cd "$BASE_DIR"
git clone https://github.com/mosp-catalysis/gcmc_kmc_Cu_CO_H cu-sim-workspace
cd cu-sim-workspace

# Compile GCMC
cd GCMC-code
gfortran -O3 GCMC.f95 -o GCMC.exe
gfortran -O3 gen-random-cov-slab.f95 -o gen-random-cov-slab.exe
cd ../KMC-code

# Compile EKMC
make
cd ..

# Compile Post-processing
cd scripts
gfortran -O3 Natom-eachlayer-inlast_one-10L.f95 -o Natom.exe
cd ..
```

## 2. Environment Setup
Inside `cu-sim-workspace`, create the specific environment for the chosen surface.

```bash
mkdir -p run_{{TARGET_SURFACE}}
cd run_{{TARGET_SURFACE}}
cp ../GCMC-code/GCMC.exe ./
cp ../GCMC-code/gen-random-cov-slab.exe ./
cp ../KMC-code/EKMC.exe ./
cp ../scripts/Natom.exe ./
cp ../initial-structures-parameters/{{TARGET_SURFACE}}/* ./
```

## 3. Universal Simulation Script Generation
Inside `run_{{TARGET_SURFACE}}`, generate and run the master script. 

```bash
cat << 'EOF' > run.sh
#!/bin/bash
set -euo pipefail

echo "=================================================="
echo "Starting simulation process (PID: $$) at $(date)"
echo "Target Surface: {{TARGET_SURFACE}}"
echo "Target Time: {{TARGET_TIME}} s"
echo "=================================================="
shopt -s nullglob
CURDIR=$PWD
TARGET_KMC_TIME={{TARGET_TIME}}

# Validate TARGET_KMC_TIME is a valid number
if ! [[ "$TARGET_KMC_TIME" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
    echo "Fatal Error: Target time '$TARGET_KMC_TIME' is not a valid number."
    exit 1
fi

REQUIRED_FILES=("GCMC.exe" "gen-random-cov-slab.exe" "EKMC.exe" "Natom.exe" "input-forGCMC" "input-forKMC1diff" "input-forKMC2aj" "ini.xyz" "Total_bulk.xyz")
for file in "${REQUIRED_FILES[@]}"; do
    if [ ! -f "$CURDIR/$file" ]; then
        echo "Fatal Error: Required file '$file' not found!"
        exit 1
    fi
done

run_gcmc() {
    local cycle=$1
    mkdir -p "${cycle}-GCMC"
    cd "${cycle}-GCMC"
    cp "${CURDIR}/GCMC.exe" ./
    cp "${CURDIR}/input-forGCMC" ./input
    cp "${CURDIR}/ini.xyz" ./
    ./GCMC.exe >> stdout
    cp "${CURDIR}/gen-random-cov-slab.exe" ./
    ./gen-random-cov-slab.exe >> stdout
    cp ./last_one_random.xyz "${CURDIR}/ini.xyz"
    cd "${CURDIR}"
}

run_kmc() {
    local cycle=$1
    local j=$2
    for step in "1diff" "2aj"; do
        mkdir -p "${cycle}-KMC-${j}-${step}"
        cd "${cycle}-KMC-${j}-${step}"
        cp "${CURDIR}/EKMC.exe" ./
        cp "${CURDIR}/ini.xyz" ./
        cp "${CURDIR}/Total_bulk.xyz" ./
        cp "${CURDIR}/input-forKMC${step}" ./input
        ./EKMC.exe >> stdout
        cp ./last_one.xyz "${CURDIR}/ini.xyz"
        cd "${CURDIR}"
    done
}

prev_cov_CO=-1000.0
prev_cov_H=-1000.0
CONVERGED_FLAG=0
# Ridiculously high safety net to prevent infinite loops if time never increases, without disrupting normal long runs
MAX_GLOBAL_CYCLES=999999  
global_cum_time=0.0
TIME_REACHED_FLAG=0
FINAL_COMPLETED_CYCLE=0

echo "Cycle Num_CO Cov_CO Num_H Cov_H" > coverage_log.txt

cycle_counter=1
while true; do
    if [ "$cycle_counter" -gt "$MAX_GLOBAL_CYCLES" ]; then
        echo "Warning: Failsafe limit of $MAX_GLOBAL_CYCLES cycles reached. Terminating."
        FINAL_COMPLETED_CYCLE=$((cycle_counter - 1))
        break
    fi

    space=$(df -m . | awk 'NR==2 {print $4}')
    if [ "$space" -lt 1000 ]; then
        echo "Fatal Error: Disk space running low (< 1000 MB). Terminating."
        FINAL_COMPLETED_CYCLE=$((cycle_counter - 1))
        break
    fi

    run_gcmc "$cycle_counter"
    
    GCMC_FILE="${cycle_counter}-GCMC/last_one.xyz"
    if [ ! -f "$GCMC_FILE" ]; then
        echo "Fatal Error: $GCMC_FILE not generated. Simulation failed."
        FINAL_COMPLETED_CYCLE=$((cycle_counter - 1))
        break
    fi

    LINE2=$(sed -n '2p' "$GCMC_FILE" || true)
    num_CO=$(echo "$LINE2" | awk '{print $4}')
    cov_CO=$(echo "$LINE2" | awk '{print $5}')
    num_H=$(echo "$LINE2" | awk '{print $6}')
    cov_H=$(echo "$LINE2" | awk '{print $7}')
    
    if [[ -z "$cov_CO" || -z "$cov_H" ]]; then
        echo "Fatal Error: Failed to parse coverage values from $GCMC_FILE."
        FINAL_COMPLETED_CYCLE=$((cycle_counter - 1))
        break
    fi

    echo "$cycle_counter $num_CO $cov_CO $num_H $cov_H" >> coverage_log.txt
    
    if [ "$cycle_counter" -gt 1 ] && [ "$CONVERGED_FLAG" -eq 0 ]; then
        is_converged=$(awk -v cCO="$cov_CO" -v pCO="$prev_cov_CO" -v cH="$cov_H" -v pH="$prev_cov_H" '
        BEGIN {
            diff_CO = (cCO - pCO < 0) ? -(cCO - pCO) : (cCO - pCO);
            diff_H = (cH - pH < 0) ? -(cH - pH) : (cH - pH);
            if (diff_CO < 0.001 && diff_H < 0.005) { print 1 } else { print 0 }
        }')
        if [ "$is_converged" -eq 1 ]; then
            echo "Convergence Reached at Cycle $cycle_counter!" >> coverage_log.txt
            CONVERGED_FLAG=1
        fi
    fi
    prev_cov_CO=$cov_CO
    prev_cov_H=$cov_H

    KMC_LIMIT=100
    [ "$CONVERGED_FLAG" -eq 1 ] && KMC_LIMIT=999999

    j=1
    while [ "$j" -le "$KMC_LIMIT" ]; do
        run_kmc "$cycle_counter" "$j"
        
        FILE_1="${cycle_counter}-KMC-${j}-1diff/last_one.xyz"
        FILE_2="${cycle_counter}-KMC-${j}-2aj/last_one.xyz"
        
        if [ ! -f "$FILE_1" ] || [ ! -f "$FILE_2" ]; then
            echo "Fatal Error: KMC files missing at cycle $cycle_counter step $j."
            FINAL_COMPLETED_CYCLE=$((cycle_counter - 1))
            break 2
        fi

        t1=$(sed -n '2p' "$FILE_1" | awk '{print $4}')
        t2=$(sed -n '2p' "$FILE_2" | awk '{print $4}')
        
        if [[ -z "$t1" || -z "$t2" ]]; then
            echo "Fatal Error: Failed to parse time values at cycle $cycle_counter step $j."
            FINAL_COMPLETED_CYCLE=$((cycle_counter - 1))
            break 2
        fi

        global_cum_time=$(awk -v g="$global_cum_time" -v a="$t1" -v b="$t2" 'BEGIN {print g+a+b}')
        
        is_time_reached=$(awk -v c="$global_cum_time" -v t="$TARGET_KMC_TIME" 'BEGIN { if (c >= t) print 1; else print 0 }')
        if [ "$is_time_reached" -eq 1 ]; then
            echo "Target KMC time reached: $global_cum_time s (>= $TARGET_KMC_TIME s) at cycle $cycle_counter step $j."
            TIME_REACHED_FLAG=1
            FINAL_COMPLETED_CYCLE=$cycle_counter
            break 2 
        fi
        j=$((j+1))
    done

    FINAL_COMPLETED_CYCLE=$cycle_counter
    cycle_counter=$((cycle_counter+1))
done

# --- POST-PROCESSING ---
echo "Starting auto post-processing at $(date)..."
> n_atoms_highest_layer_tot.dat

for (( c=1; c<=FINAL_COMPLETED_CYCLE; c++ )); do
    k=1
    while true; do
        f1="${c}-KMC-${k}-1diff"
        f2="${c}-KMC-${k}-2aj"
        if [ ! -d "$f1" ]; then break; fi

        cd "$f1" || exit 1
        cp "${CURDIR}/Natom.exe" ./ && ./Natom.exe > /dev/null
        [ -f ./n_atoms_highest_layer.dat ] && echo "$f1  $(cat ./n_atoms_highest_layer.dat)" >> "${CURDIR}/n_atoms_highest_layer_tot.dat"
        cd "${CURDIR}"

        if [ -d "$f2" ]; then
            cd "$f2" || exit 1
            cp "${CURDIR}/Natom.exe" ./ && ./Natom.exe > /dev/null
            [ -f ./n_atoms_highest_layer.dat ] && echo "$f2  $(cat ./n_atoms_highest_layer.dat)" >> "${CURDIR}/n_atoms_highest_layer_tot.dat"
            cd "${CURDIR}"
        fi
        k=$((k+1))
    done
done

# --- PURE PYTHON PROCESSING ---
cat << 'PYEOF' > process_data.py
import sys
import os
import csv

root_path = sys.argv[1]
input_file = sys.argv[2]
output_file = sys.argv[3]
file_path = os.path.join(root_path, input_file)

# First pass: read all valid rows and determine the global maximum width
raw_rows = []
cum_sums = [0.0, 0.0, 0.0]
max_cols = 4  # Folder name + 3 cumulative columns at minimum

if os.path.exists(file_path):
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 4:
                continue
            try:
                c1 = float(parts[1])
                c2 = float(parts[2])
                c3 = float(parts[3])
            except ValueError:
                continue

            cum_sums[0] += c1
            cum_sums[1] += c2
            cum_sums[2] += c3
            max_cols = max(max_cols, len(parts))

            # Store the current row in the raw form first
            raw_rows.append([parts[0], cum_sums[0], cum_sums[1], cum_sums[2]] + parts[4:])

# Second pass: build a header and pad every row to the same maximum width
if raw_rows:
    header = ["Folder_Name", "Cum_Total_Steps", "Cum_Migration_Steps", "Cum_Time_s"]

    extra_cols = max_cols - 4
    for i in range(extra_cols):
        pair_idx = i // 2 + 1
        if i % 2 == 0:
            if pair_idx == 1:
                header.append("Z_Layer_1(Top)")
            else:
                header.append("Z_Layer_{}".format(pair_idx))
        else:
            header.append("Atom_Count_{}".format(pair_idx))

    # Pad header and rows to the same width
    if len(header) < max_cols:
        for i in range(len(header), max_cols):
            header.append("Column_{}".format(i + 1))

    padded_rows = []
    for row in raw_rows:
        if len(row) < max_cols:
            row = row + [""] * (max_cols - len(row))
        elif len(row) > max_cols:
            row = row[:max_cols]
        padded_rows.append(row)

    with open(os.path.join(root_path, output_file), 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(padded_rows)
PYEOF

python process_data.py "${CURDIR}" n_atoms_highest_layer_tot.dat KMC_accumulated_results.csv
echo "Post-processing complete. Results saved to KMC_accumulated_results.csv"

echo "=================================================="
echo "Simulation process (PID: $$) done at $(date)"
echo "=================================================="
EOF

chmod +x run.sh
nohup ./run.sh > simulation_master.log 2>&1 &
```

Inform the user that the job has been launched in the background successfully. Tell them the exact directory (`run_{{TARGET_SURFACE}}`) and that they can monitor `simulation_master.log`.