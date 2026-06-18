#!/bin/bash
#SBATCH -p cp4
#SBATCH -o testall.out
#SBATCH -J pyscf-test
#SBATCH -n 56

start=$(date +%s)
echo "Job started on $(date)"

# ================= Task =================
python DyCp2.py -i ./
# =======================================

# Execute the following commands only after the task has finished
echo "Task finished, cleaning .h5 files..."

# Get the directory where this shell script is located
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

# Delete only the .h5 files in that directory, without recursion
rm -f "$SCRIPT_DIR"/*.h5

echo "Cleanup done."

end=$(date +%s)
echo "Job finished on $(date)"
echo "Elapsed time: $((end - start)) seconds"
