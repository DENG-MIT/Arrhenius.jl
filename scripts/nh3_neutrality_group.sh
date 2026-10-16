#!/usr/bin/env bash
set -Eeuo pipefail
[[ $# == 4 ]] || { echo "usage: nh3_neutrality_group.sh <Arrhenius-project> <ids.txt> <output-dir> <nh3-parity-root>" >&2; exit 2; }
project=$(realpath "$1"); ids=$(realpath "$2"); out=$(realpath -m "$3"); root=$(realpath "$4")
[[ -f $project/Project.toml && -f $ids && -f $root/validation_set.txt ]] || { echo "missing project, id list, or validation set" >&2; exit 2; }
source "$root/env.sh"
unset NH3_DISCRETIZATION
mkdir -p "$out"
manifest="$out/manifest.csv"
"$PY" "$root/h/nh3_set_manifest.py" "$root" "$manifest" >/dev/null
scratch=$(mktemp -d "${TMPDIR:-/tmp}/yamlnative-neutrality.XXXXXX")
fifo=$scratch/request.fifo; mkfifo "$fifo"
server_pid=
cleanup() {
    local rc=$?
    if [[ -n ${server_pid:-} ]] && kill -0 "$server_pid" 2>/dev/null; then
        printf 'QUIT\n' >&9 2>/dev/null || true
        timeout -k 10 20 tail --pid="$server_pid" -f /dev/null 2>/dev/null || kill "$server_pid" 2>/dev/null || true
    fi
    rm -rf "$scratch"
    exit "$rc"
}
trap cleanup EXIT INT TERM
cd "$root"
"$JL" --project="$project" "$root/h/nh3_native_server.jl" < "$fifo" > "$out/native.log" 2>&1 &
server_pid=$!
exec 9>"$fifo"
count=0
while IFS= read -r id || [[ -n $id ]]; do
    [[ -n $id ]] || continue
    spec=$(awk -F, -v id="$id" '$1 == id {print $2}' "$manifest")
    [[ -n $spec ]] || { echo "unknown validation id: $id" >&2; exit 2; }
    marker=$scratch/$id.done; output=$out/$id.npz
    printf '%s\t%s\t0\t0\t%s\n' "$root/$spec" "$output" "$marker" >&9
    status=$(timeout -k 10 900 "$root/h/wait_marker.sh" "$marker") || {
        tail -80 "$out/native.log" >&2; echo "case $id did not finish within 900 seconds" >&2; exit 124;
    }
    [[ $status == DONE ]] || { tail -80 "$out/native.log" >&2; echo "case $id failed: $status" >&2; exit 1; }
    rm -f "$marker"
    count=$((count+1)); echo "completed $id ($count)"
done < "$ids"
printf 'QUIT\n' >&9
exec 9>&-
wait "$server_pid"; server_pid=
echo "group complete: $count cases project=$project output=$out"