#!/bin/bash
# Launch every experiment that backs a manuscript FIGURE.
#
# EXCLUDED: the two Jupyter-notebook studies -- appendix figures A1-A6
# (exp_simplified_methods.ipynb) and A13-A14 (exp_richardson.ipynb) -- which are
# run from the notebooks themselves, not from a submit script.
#
# Each line below runs one submit_experiments_*.sh with a specific CONF. Those
# scripts submit one `sbatch` job per configuration and random seed, so this is
# meant to be run from the head node of a Slurm cluster. Jobs whose output .txt
# already exists are skipped, so it is safe to re-run.
#
# Usage:
#   ./submit_paper_experiments.sh            # submit everything
#   ./submit_paper_experiments.sh --dry-run  # print what would run, submit nothing
#   ./submit_paper_experiments.sh --list     # show the figure -> (script, CONF) map
#
# The CIFAR-10 (A31) and BigEarthNet (Fig 6) rows need the real data downloaded
# first -- see ../data/README.md (BigEarthNet is > 30 GB).

set -u
cd "$(dirname "$0")"

# figure(s) | launcher | CONF
RUNS=(
  "Figures 1, A16, A17|submit_experiments_1.sh|1"
  "Figures 2, A18|submit_experiments_1.sh|2"
  "Figures 3, A19-A23|submit_experiments_1.sh|3"
  "Figure 4|submit_experiments_1.sh|4"
  "Figure 5|submit_experiments_4.sh|401"
  "Figure 6 (needs data/bigearthnet-full)|submit_experiments_6_bigearthnet.sh|601"
  "Figures A7, A9, A11|submit_experiments_2_comparisons.sh|203"
  "Figures A8, A10, A12|submit_experiments_2_comparisons.sh|204"
  "Figure A24|submit_experiments_2_comparisons.sh|201"
  "Figures A25, A26|submit_experiments_4.sh|403"
  "Figure A27|submit_experiments_4.sh|406"
  "Figure A28|submit_experiments_3_T_estimation.sh|301"
  "Figure A29|submit_experiments_3_T_estimation.sh|304"
  "Figure A30|submit_experiments_3_T_estimation.sh|305"
  "Figure A31 (needs data/cifar10, data/cifar10h)|submit_experiments_5_cifar10.sh|503"
)

case "${1:-}" in
  --list)
    printf '%-42s  %-38s  %s\n' "FIGURE(S)" "LAUNCHER" "CONF"
    for entry in "${RUNS[@]}"; do
      IFS='|' read -r figs script conf <<< "$entry"
      printf '%-42s  %-38s  %s\n' "$figs" "$script" "$conf"
    done
    exit 0 ;;
  --dry-run) DRYRUN=1 ;;
  "")        DRYRUN=0 ;;
  *) echo "unknown option: $1  (use --dry-run or --list)" >&2; exit 1 ;;
esac

failed=()
for entry in "${RUNS[@]}"; do
  IFS='|' read -r figs script conf <<< "$entry"
  echo
  echo "=================================================================="
  echo ">> $figs   ->   CONF=$conf  $script"
  echo "=================================================================="
  if [ "$DRYRUN" -eq 1 ]; then
    echo "   (dry run)  CONF=$conf ./$script"
  else
    if ! CONF="$conf" bash "./$script"; then
      echo "   !! $script (CONF=$conf) exited non-zero" >&2
      failed+=("$script CONF=$conf")
    fi
  fi
done

echo
if [ "${#failed[@]}" -ne 0 ]; then
  echo "Finished with errors in:"
  printf '  %s\n' "${failed[@]}"
  exit 1
fi
echo "All figure experiments submitted. Re-run any time; completed jobs are skipped."
echo "Then collect results into experiments/results/ and build figures with make_plots.R"
echo "(see ../metadata/workflow.md)."
