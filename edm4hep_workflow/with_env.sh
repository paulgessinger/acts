#!/bin/bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}"   )" &> /dev/null && pwd   )

source ~/spack/share/spack/setup-env.sh
spack env activate ~/dev/ci-dependencies/

source "$SCRIPT_DIR/../build/this_acts_withdeps.sh"

source "$SCRIPT_DIR/.venv/bin/activate"

exec "$@"
