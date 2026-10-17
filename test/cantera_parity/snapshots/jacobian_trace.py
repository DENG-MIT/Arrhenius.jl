"""Parse Cantera level-2 Newton logs with per-attempt and cumulative Jacobian counts."""
import re


def latest_stat(stats, key):
    """Read Cantera's cumulative Jacobian counter snapshot at the last save."""
    values = stats.get(key, [])
    return int(values[-1]) if values else 0


def cantera_attempt_trace(path):
    rows = []
    current = None
    last_jac = 0
    pending_refreshes = 0
    with open(path, encoding="utf-8", errors="replace") as stream:
        for line in stream:
            match = re.search(r"Timestep \(\d+\) dt=\s*([0-9.eE+-]+)", line)
            if match:
                if current is not None:
                    rows.append(current)
                current = {"dt": float(match.group(1)), "damping": [],
                           "accepted": None, "reset_after": False,
                           "jacobians": 0, "jacobians_cumulative": last_jac,
                           "jacobian_age_start": None}
                continue
            if current is not None and ("Re-evaluating Jacobian" in line or
                                        "Maximum Jacobian age reached" in line):
                current["jacobians"] += 1
                pending_refreshes += 1
                if (current["jacobian_age_start"] is None and
                        "Maximum Jacobian age reached" in line):
                    current["jacobian_age_start"] = 0
            columns = line.split()
            if (len(columns) == 8 and columns[0].isdigit() and
                    re.fullmatch(r"[0-9.eE+-]+", columns[1]) and "/" in columns[-1]):
                njac = int(columns[6])
                age = int(columns[7].split("/", 1)[0])
                delta = max(0, njac - last_jac)
                if current is not None:
                    newly_observed = max(0, delta - pending_refreshes)
                    current["jacobians"] += newly_observed
                    current["jacobians_cumulative"] = njac
                    current["damping"].append(float(columns[1]))
                    if current["jacobian_age_start"] is None:
                        current["jacobian_age_start"] = 0 if newly_observed else age - 1
                pending_refreshes = max(0, pending_refreshes - delta)
                last_jac = njac
            if current is None:
                continue
            if "Resetting negative species" in line:
                current["reset_after"] = True
            if re.search(r"Timestep \(\d+\) failed", line):
                current["accepted"] = False
            elif re.search(r"Timestep \(\d+\) succeeded", line):
                current["accepted"] = True
                rows.append(current)
                current = None
    if current is not None:
        rows.append(current)
    return rows
