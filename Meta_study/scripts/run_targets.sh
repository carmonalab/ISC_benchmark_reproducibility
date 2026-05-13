#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."
Rscript -e "targets::tar_make()"
