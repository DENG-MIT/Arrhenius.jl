"""Host checks shared by the paired calculation benchmarks."""
import hashlib
import ctypes
import os
from pathlib import Path
import platform
import subprocess


def cpu_brand():
    if platform.system() == "Darwin":
        return subprocess.check_output(
            ["sysctl", "-n", "machdep.cpu.brand_string"], text=True).strip()
    if platform.system() == "Linux":
        for line in Path("/proc/cpuinfo").read_text().splitlines():
            if line.startswith("model name"):
                return line.split(":", 1)[1].strip()
    return platform.processor() or "unknown"


def host_metadata():
    return {"cpu": cpu_brand(), "platform": platform.platform(),
            "kernel_release": platform.release(), "system": platform.system(),
            "machine": platform.machine(), "logical_cpus": os.cpu_count()}


def matches_target(host, target):
    if target == "wsl":
        return (host.get("system") == "Linux"
                and "microsoft" in host.get("kernel_release", "").lower())
    if target == "apple-m4":
        return host.get("system") == "Darwin" and "Apple M4" in host.get("cpu", "")
    raise ValueError(f"unknown benchmark target: {target}")


def loaded_library_paths():
    """Return mapped library paths on the two supported benchmark platforms."""
    paths = set()
    maps = Path("/proc/self/maps")
    if maps.is_file():
        for line in maps.read_text().splitlines():
            path = line.split()[-1]
            if path.startswith("/"):
                paths.add(Path(path).resolve())
    elif platform.system() == "Darwin":
        library = ctypes.CDLL(None)
        count = library._dyld_image_count
        count.argtypes, count.restype = [], ctypes.c_uint32
        name = library._dyld_get_image_name
        name.argtypes, name.restype = [ctypes.c_uint32], ctypes.c_char_p
        for index in range(count()):
            path = name(index)
            if path:
                paths.add(Path(os.fsdecode(path)).resolve())
    return sorted(paths)


def cantera_library_hashes(module_path):
    """Hash the actual mapped Cantera libraries; module_path is retained for callers."""
    paths = {p for p in loaded_library_paths() if "libcantera" in p.name and p.is_file()}
    return {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(paths)}


def verify_numerical_threads(*, set_accelerate=False):
    """Check loaded BLAS/OpenMP runtimes directly, without optional packages."""
    pools = []
    for path in loaded_library_paths():
        name = path.name.lower()
        if "openblas" in name:
            symbols = ("openblas_get_num_threads", "openblas_get_num_threads64_",
                       "scipy_openblas_get_num_threads", "scipy_openblas_get_num_threads64_")
        elif "mkl_rt" in name:
            symbols = ("MKL_Get_Max_Threads", "mkl_get_max_threads")
        elif any(part in name for part in ("libgomp", "libiomp", "libomp")):
            symbols = ("omp_get_max_threads",)
        else:
            continue
        library = ctypes.CDLL(str(path))
        for symbol in symbols:
            query = getattr(library, symbol, None)
            if query is not None:
                query.argtypes, query.restype = [], ctypes.c_int
                pools.append({"file": path.name, "query": symbol, "threads": int(query())})
                break
        else:
            raise RuntimeError(f"cannot verify loaded numerical runtime {path}")
    if not pools or any(pool["threads"] != 1 for pool in pools):
        raise RuntimeError(f"expected verified single-thread numerical runtimes: {pools}")
    result = {"threadpools": pools}
    if platform.system() == "Darwin":
        library = ctypes.CDLL("/System/Library/Frameworks/Accelerate.framework/Accelerate")
        getter = library.BLASGetThreading
        getter.argtypes, getter.restype = [], ctypes.c_uint
        if set_accelerate:
            setter = library.BLASSetThreading
            setter.argtypes, setter.restype = [ctypes.c_uint], ctypes.c_int
            if setter(1) != 0:
                raise RuntimeError("Accelerate single-thread setting failed")
        result["accelerate_threading_mode"] = int(getter())
        if result["accelerate_threading_mode"] != 1:
            raise RuntimeError("Accelerate is not set to single-thread mode")
    return result
