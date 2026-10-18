#!/usr/bin/env bash
set -Eeuo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
IMAGE="${PROJECT_DIR}/containers/images/srebf2-phase12_0.1.sif"
VALIDATION_DIR="${PROJECT_DIR}/singularity_validation/phase12"
REPORT_DIR="${PROJECT_DIR}/reports/workflow"

mkdir -p "${REPORT_DIR}"

echo "=== SREBF2 Phase 12 Singularity Validation ==="
echo "Date: $(date -Is)"
echo "Git commit: $(git -C "${PROJECT_DIR}" rev-parse --short HEAD)"

if [[ ! -f "${IMAGE}" ]]; then
    echo "ERROR: Singularity image not found: ${IMAGE}"
    exit 1
fi

if command -v singularity >/dev/null 2>&1; then
    RUNTIME=singularity
elif command -v apptainer >/dev/null 2>&1; then
    RUNTIME=apptainer
else
    echo "ERROR: Singularity/Apptainer is not installed."
    exit 1
fi

"${RUNTIME}" --version
sha256sum "${IMAGE}"

echo "=== Preparing isolated workspace ==="

if [[ -e "${VALIDATION_DIR}" ]]; then
    echo "ERROR: Validation workspace already exists:"
    echo "${VALIDATION_DIR}"
    echo "Move or archive it before starting a new validation."
    exit 1
fi

mkdir -p \
    "${VALIDATION_DIR}/scripts/phase12" \
    "${VALIDATION_DIR}/results/phase12" \
    "${VALIDATION_DIR}/reports/phase12"

cp "${PROJECT_DIR}/scripts/phase12/phase12_robustness.py" \
   "${VALIDATION_DIR}/scripts/phase12/"

cp "${PROJECT_DIR}/results/phase12/SREBF2_HELA_REGULATORY_SIGNAL.tsv" \
   "${VALIDATION_DIR}/results/phase12/"

echo "=== Checking input integrity ==="

cmp \
  "${PROJECT_DIR}/results/phase12/SREBF2_HELA_REGULATORY_SIGNAL.tsv" \
  "${VALIDATION_DIR}/results/phase12/SREBF2_HELA_REGULATORY_SIGNAL.tsv"

echo "Input matrix: IDENTICAL"

echo "=== Checking scientific dependencies ==="

"${RUNTIME}" exec --cleanenv "${IMAGE}" python -c '
import numpy, pandas, scipy, statsmodels, patsy, pyBigWig

print("numpy:", numpy.__version__)
print("pandas:", pandas.__version__)
print("scipy:", scipy.__version__)
print("statsmodels:", statsmodels.__version__)
print("patsy:", patsy.__version__)
print("pyBigWig:", pyBigWig.__version__)
print("Scientific imports: PASS")
'

echo "=== Running Phase 12 robustness analysis ==="

"${RUNTIME}" exec \
    --cleanenv \
    --bind "${VALIDATION_DIR}:/workspace" \
    --pwd /workspace \
    "${IMAGE}" \
    python scripts/phase12/phase12_robustness.py \
    > "${VALIDATION_DIR}/reports/phase12/phase12_robustness.log" 2>&1

echo "=== Comparing statistical outputs ==="

PROJECT_DIR="${PROJECT_DIR}" \
VALIDATION_DIR="${VALIDATION_DIR}" \
python3 - <<'PY'
import os
from pathlib import Path

import numpy as np
import pandas as pd

project = Path(os.environ["PROJECT_DIR"])
validation = Path(os.environ["VALIDATION_DIR"])

native = project / "results/phase12"
singularity = validation / "results/phase12"

files = [
    "PHASE12_ROBUSTNESS_MODELS.tsv",
    "PHASE12_VIF.tsv",
]

records = []

for name in files:
    a = pd.read_csv(native / name, sep="\t")
    b = pd.read_csv(singularity / name, sep="\t")

    structure_ok = (
        a.shape == b.shape
        and a.columns.equals(b.columns)
    )

    equivalent = structure_ok
    max_difference = 0.0

    if structure_ok:
        for column in a.columns:
            if (
                pd.api.types.is_numeric_dtype(a[column])
                and pd.api.types.is_numeric_dtype(b[column])
            ):
                x = a[column].to_numpy(dtype=float)
                y = b[column].to_numpy(dtype=float)

                equivalent &= bool(
                    np.allclose(
                        x,
                        y,
                        rtol=1e-8,
                        atol=1e-10,
                        equal_nan=True,
                    )
                )

                finite = np.isfinite(x) & np.isfinite(y)

                if finite.any():
                    max_difference = max(
                        max_difference,
                        float(np.max(np.abs(x[finite] - y[finite])))
                    )
            else:
                equivalent &= bool(a[column].equals(b[column]))

    records.append({
        "file": name,
        "structure_identical": structure_ok,
        "numerically_equivalent": equivalent,
        "maximum_absolute_difference": max_difference,
        "rtol": 1e-8,
        "atol": 1e-10,
    })

result = pd.DataFrame(records)

output = (
    project /
    "reports/workflow/phase12_singularity_reproducibility.tsv"
)

result.to_csv(output, sep="\t", index=False)

print(result.to_string(index=False))
print(f"\nSaved: {output}")

if not result["numerically_equivalent"].all():
    raise SystemExit("Singularity reproducibility: FAIL")

print("\nSingularity reproducibility: PASS")
PY

echo "=== PHASE 12 SINGULARITY VALIDATION: PASS ==="
