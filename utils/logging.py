"""Opt-in logging utilities: tee stdout/stderr to a timestamped log file."""

import os
import sys
import time


class _Tee:
    """Write to both a file and the original stream."""
    def __init__(self, stream, logfile):
        self._stream = stream
        self._log = logfile

    def write(self, msg):
        self._stream.write(msg)
        self._log.write(msg)
        self._log.flush()

    def flush(self):
        self._stream.flush()
        self._log.flush()


def init_logging(log_dir=None):
    """Set up stdout/stderr tee to a timestamped log file.

    Parameters
    ----------
    log_dir : str or None
        Directory for log files. Defaults to ``logs/`` under the caller's
        script directory.

    Returns
    -------
    str
        Path to the log file.
    """
    if log_dir is None:
        log_dir = os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])), "logs")
    os.makedirs(log_dir, exist_ok=True)

    ts = time.strftime("%Y%m%d_%H%M%S")
    log_path = os.path.join(log_dir, f"run_{ts}.log")

    fh = open(log_path, "w")
    fh.write(f"# {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
    fh.write(f"# Args: {sys.argv}\n\n")

    sys.stdout = _Tee(sys.__stdout__, fh)
    sys.stderr = _Tee(sys.__stderr__, fh)
    print(f"Logging to: {log_path}")
    return log_path
