import sys

"""
Version has unique source in pyproject.toml.
importlib fetches version from distribution metadata files
(in dist-info or egg-info dirs).
From Python 3.8, importlib_metadata is in standard library as importlib.metadata.
"""
if sys.version_info >= (3, 8):
    from importlib import metadata
else:
    import importlib_metadata as metadata

try:
    __version__: str = metadata.version("pling")
except metadata.PackageNotFoundError:
    __version__ = "0+unknown"


def main():
    from .run_pling import main as run_main

    return run_main()


__all__ = ["main"]
