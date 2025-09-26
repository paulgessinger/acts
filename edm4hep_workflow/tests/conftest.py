"""Pytest configuration and shared fixtures."""

import pytest
from pathlib import Path


@pytest.fixture(scope="session")
def project_root():
    """Return the project root directory."""
    return Path(__file__).parent.parent


@pytest.fixture
def sample_data_dir(project_root):
    """Return path to sample data directory."""
    return project_root / "tests" / "data"
