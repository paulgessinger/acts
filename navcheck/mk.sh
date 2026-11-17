#!/usr/bin/env bash
# Usage: mk.sh <acts worktree> <output binary> [extra compiler flags...]
# Compiles navcheck.cpp against the build tree of the given worktree, reusing
# the flags of one of its JSON plugin translation units.
set -euo pipefail
tree=$1
out=$2
shift 2
here=$(cd "$(dirname "$0")" && pwd)
flags=$(python3 -I - "$tree" <<'EOF'
import json, shlex, sys
tree = sys.argv[1]
cc = json.load(open(f"{tree}/build/compile_commands.json"))
entry = next(e for e in cc if e["file"].endswith("TrackingGeometryJsonConverter.cpp"))
args = shlex.split(entry["command"])
keep = []
skip = False
for a in args[1:]:
    if skip:
        skip = False
        continue
    if a in ("-o", "-c", "-MF", "-MT"):
        skip = True
        continue
    if a in ("-MD",) or a.endswith(".cpp") or a.endswith("ccache"):
        continue
    if a.startswith("-Werror"):
        continue
    keep.append(a)
print(" ".join(shlex.quote(a) for a in keep))
EOF
)
cd "$tree"
eval ./build.sh c++ $flags -I"$tree/Plugins/Json/include" "$@" \
  "$here/navcheck.cpp" -o "$out" \
  -L"$tree/build/lib" -lActsCore -lActsPluginJson \
  -Wl,-rpath,"$tree/build/lib"
