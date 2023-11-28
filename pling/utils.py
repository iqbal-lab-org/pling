from pathlib import Path


def get_pling_root_dir() -> Path:
    return Path(__file__).parent.parent
