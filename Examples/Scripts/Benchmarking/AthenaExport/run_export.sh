#!/usr/bin/env bash
# Run from a clean shell. Build the exporter first; see README.md.
set -euo pipefail
workdir=${1:?Usage: run_export.sh OUTPUT_DIRECTORY [EVENTS]}
events=${2:-2}
if [[ -e "${workdir}/tracks" || -e "${workdir}/AOD.pool.root" ]]; then
  echo "Use a fresh output directory to avoid mixing event samples: ${workdir}" >&2
  exit 1
fi
athena_release=${AMVF_ATHENA_RELEASE:-/cvmfs/atlas-nightlies.cern.ch/repo/sw/main_Athena_x86_64-el9-gcc15-opt/2026-09-14T2100/Athena/25.0.73/InstallArea/x86_64-el9-gcc15-opt}
export_build=${AMVF_EXPORT_BUILD:-/tmp/amvf-athena-export}
# Release setup scripts are not compatible with nounset.
set +u
source /cvmfs/sft.cern.ch/lcg/releases/gcc/15.2.0/x86_64-el9/setup.sh
source "${athena_release}/setup.sh"
source "${export_build}/x86_64-el9-gcc15-opt/setup.sh"
set -u
# Standard fallback from ATLAS asetup's epilog; preserve site configuration.
export FRONTIER_SERVER="${FRONTIER_SERVER:-(serverurl=http://atlasfrontier-ai.cern.ch:8000/atlr)(proxyurl=http://v4f.hl-lhc.net:6082)}"
export ATLAS_POOLCOND_PATH="${ATLAS_POOLCOND_PATH:-/cvmfs/atlas-condb.cern.ch/repo/conditions}"
mkdir -p "${workdir}"
cd "${workdir}"
export AMVF_EXPORT_DIR="${PWD}/tracks"
export ATHENA_CORE_NUMBER=1
rdo=$(python -c 'from AthenaConfiguration.TestDefaults import defaultTestFiles; print(defaultTestFiles.RDO_RUN4[0])')
conditions=$(python -c 'from AthenaConfiguration.TestDefaults import defaultConditionsTags; print(defaultConditionsTags.RUN4_MC)')
python - "$rdo" "$conditions" "$athena_release" <<'PY'
import json, sys
with open('export-provenance.json', 'w') as f:
    json.dump(dict(input_rdo=sys.argv[1], conditions=sys.argv[2], athena_release=sys.argv[3],
                   preinclude=['InDetConfig.ConfigurationHelpers.OnlyTrackingRecoPreInclude',
                               'ActsConfig.ActsCIFlags.actsProductionFlags']), f, indent=2)
PY
Reco_tf.py \
  --conditionsTag "default:${conditions}" \
  --preInclude "InDetConfig.ConfigurationHelpers.OnlyTrackingRecoPreInclude,ActsConfig.ActsCIFlags.actsProductionFlags" \
  --postInclude "RAWtoALL:AmvfExport.ExportConfig.addExport" \
  --inputRDOFile "${rdo}" --outputAODFile AOD.pool.root --maxEvents "${events}"
