#!/usr/bin/env bash
set -euo pipefail
repo=$1; manifest=$2; out=$3; batches=$4; per_batch=$5
cd "${NH3_HARNESS_DIR:?set NH3_HARNESS_DIR to the validation harness directory}"
source ./env.sh
export LC_ALL=C OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p "$out"
scratch=$(mktemp -d "${TMPDIR:-/tmp}/kernel-neutrality.XXXXXX")
mkfifo "$scratch/request.fifo"
server_pid=
cleanup() {
    local rc=$?
    if [[ -n "$server_pid" ]] && kill -0 "$server_pid" 2>/dev/null; then
        printf 'QUIT\n' >&9 2>/dev/null || true
        timeout 30 tail --pid="$server_pid" -f /dev/null 2>/dev/null || kill "$server_pid" 2>/dev/null || true
    fi
    rm -rf "$scratch"
    exit "$rc"
}
trap cleanup EXIT INT TERM
"$JL" --project="$repo" h/nh3_native_server.jl < "$scratch/request.fifo" > "$scratch/native.log" 2>&1 &
server_pid=$!
exec 9>"$scratch/request.fifo"
while IFS=, read -r id spec; do
    [[ -n "$id" && -n "$spec" ]] || continue
    done_marker="$scratch/$id.done"
    printf '%s\t%s\t%s\t%s\t%s\n' "$spec" "$out/$id.npz" "$batches" "$per_batch" "$done_marker" >&9
    deadline=$((SECONDS + 900))
    while [[ ! -f "$done_marker" ]]; do
        (( SECONDS < deadline )) || { echo "timeout waiting for $id" >&2; cat "$scratch/native.log" >&2; exit 124; }
        sleep 1
    done
    status=$(cat "$done_marker")
    rm -f "$done_marker"
    [[ "$status" == DONE ]] || { echo "$id: $status" >&2; cat "$scratch/native.log" >&2; exit 1; }
done < "$manifest"
printf 'QUIT\n' >&9
exec 9>&-
wait "$server_pid"
server_pid=
echo "completed $(wc -l < "$manifest") requests for $repo"
