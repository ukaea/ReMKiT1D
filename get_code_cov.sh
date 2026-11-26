#!/bin/bash

# This file is part of ReMKiT1D.
#
# ReMKiT1D is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#
# ReMKiT1D is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with ReMKiT1D. If not, see <https://www.gnu.org/licenses/>.
#
# Copyright 2023 United Kingdom Atomic Energy Authority (stefan.mijin@ukaea.uk)
#
set -euo pipefail

LOCATION="$1"
ROOT_PATH="$2"
GCOV_BIN="${3:-gcov}"

BUILD_DIR="$(pwd)"
mkdir -p "${LOCATION}"

find "${ROOT_PATH}" \( -name '*.f90' -o -name '*.F90' \) -exec cp -t "${LOCATION}" {} +
find "${BUILD_DIR}" -path "${LOCATION}" -prune -o -name '*.gcno' -exec cp -t "${LOCATION}" {} +
find "${BUILD_DIR}" -path "${LOCATION}" -prune -o -name '*.gcda' -exec cp -t "${LOCATION}" {} +

cd "${LOCATION}"
find . -maxdepth 1 -name '*.gcda' -exec "${GCOV_BIN}" -pb {} +