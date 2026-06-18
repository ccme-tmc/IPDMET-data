#!/bin/bash
#SBATCH -p cp6
#SBATCH -o testall.out
#SBATCH -J pyscf-test
#SBATCH -n 56

start=$(date +%s)
echo "Job started on $(date)"

echo "SLURM_SUBMIT_DIR = ${SLURM_SUBMIT_DIR:-<unset>}"
echo "PWD(before task) = $(pwd)"

python CePhMe.py -i ./
rc=$?

echo "Task finished with exit code: $rc"
echo "PWD(after task)  = $(pwd)"

echo "Task finished, cleaning .h5 files..."

# Use the submission directory as the target; fall back to the current directory if unset
TARGET_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"
echo "Cleaning directory: $TARGET_DIR"

echo "Files to delete:"
find "$TARGET_DIR" -maxdepth 1 -type f -name "*.h5" -print

find "$TARGET_DIR" -maxdepth 1 -type f -name "*.h5" -delete

echo "Cleanup done."

end=$(date +%s)
echo "Job finished on $(date)"
echo "Elapsed time: $((end - start)) seconds"
exit $rc





