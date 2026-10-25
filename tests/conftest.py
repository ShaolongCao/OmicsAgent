import os

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("NUMBA_NUM_THREADS", "2")

import pytest

from omicsagent.demo import make_demo


@pytest.fixture
def demo_dir(tmp_path):
    root = tmp_path / "data"
    make_demo(root)
    return root
