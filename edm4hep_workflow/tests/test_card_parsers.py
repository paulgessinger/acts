"""Tests for MadGraph card parsing and customization functions."""

import pytest
from pathlib import Path
import shutil
from src.colliderml.cli.madgraph import (
    parse_madgraph_card,
    update_madgraph_card,
    parse_pythia8_card,
    update_pythia8_card,
)


@pytest.fixture
def test_data_dir():
    """Return path to test data directory."""
    return Path(__file__).parent / "data"


# Using pytest's built-in tmp_path fixture instead of custom temp_dir


class TestMadGraphCardParser:
    """Tests for MadGraph run_card.dat parsing."""

    def test_parse_run_card(self, test_data_dir):
        """Test parsing of real run_card.dat file."""
        run_card_path = test_data_dir / "run_card.dat"
        assert run_card_path.exists(), "run_card.dat test data not found"

        parsed = parse_madgraph_card(run_card_path)

        # Check that we parsed some expected keys
        assert isinstance(parsed, dict)
        assert len(parsed) > 0

        # Check specific keys that should be in a run card
        expected_keys = ["nevents", "iseed", "ebeam1", "ebeam2"]
        for key in expected_keys:
            assert key in parsed, f"Expected key '{key}' not found in parsed data"

        # Check that values are strings
        for key, value in parsed.items():
            assert isinstance(value, str), f"Value for key '{key}' should be string"

    def test_parse_run_card_from_string(self):
        """Test parsing run card from string content."""
        test_content = """
#*********************************************************************
# Run card
#*********************************************************************
  tag_1     = run_tag ! name of the run 
  10000 = nevents ! Number of unweighted events requested 
  0   = iseed   ! rnd seed (0=assigned automatically=default))
     1        = lpp1    ! beam 1 type 
     6500.0     = ebeam1  ! beam 1 total energy in GeV
"""
        parsed = parse_madgraph_card(test_content)

        expected = {
            "run_tag": "tag_1",
            "nevents": "10000",
            "iseed": "0",
            "lpp1": "1",
            "ebeam1": "6500.0",
        }

        for key, expected_val in expected.items():
            assert key in parsed
            assert parsed[key] == expected_val


class TestMadGraphCardUpdate:
    """Tests for MadGraph run_card.dat updating."""

    def test_update_run_card(self, test_data_dir, tmp_path):
        """Test updating a run card with new values."""
        # Copy original to temp directory
        original_path = test_data_dir / "run_card.dat"
        test_path = tmp_path / "run_card.dat"
        shutil.copy(original_path, test_path)

        # Parse original values
        original_parsed = parse_madgraph_card(test_path)
        original_nevents = original_parsed.get("nevents")

        # Update with new values
        updates = {"nevents": "5000", "iseed": "42", "new_param": "123"}

        update_madgraph_card(test_path, updates, tmp_path)

        # Parse updated values
        updated_parsed = parse_madgraph_card(test_path)

        # Check updates were applied
        assert updated_parsed["nevents"] == "5000"
        assert updated_parsed["iseed"] == "42"

        # Check if new parameters were added (may depend on implementation)
        if "new_param" in updated_parsed:
            assert updated_parsed["new_param"] == "123"

        # Check that other values weren't changed
        for key, value in original_parsed.items():
            if key not in updates:
                assert updated_parsed[key] == value


class TestPythia8CardParser:
    """Tests for Pythia8 card parsing."""

    def test_parse_pythia8_card(self, test_data_dir):
        """Test parsing of real pythia8_card_default.dat file."""
        pythia8_card_path = test_data_dir / "pythia8_card_default.dat"
        assert (
            pythia8_card_path.exists()
        ), "pythia8_card_default.dat test data not found"

        parsed = parse_pythia8_card(pythia8_card_path)

        # Check that we parsed some data
        assert isinstance(parsed, dict)
        assert len(parsed) > 0

        # Check that values are strings
        for key, value in parsed.items():
            assert isinstance(key, str)
            assert isinstance(value, str)

    def test_parse_pythia8_card_from_string(self):
        """Test parsing pythia8 card from string content."""
        test_content = """
# Pythia8 configuration
Main:numberOfEvents = 1000  # number of events
Beams:idA = 2212  # proton
Beams:idB = 2212  # proton  
Beams:eCM = 13000  # center-of-mass energy
# This is a comment line
PartonLevel:ISR = on  # initial state radiation
"""
        parsed = parse_pythia8_card(test_content)

        expected = {
            "Main:numberOfEvents": "1000",
            "Beams:idA": "2212",
            "Beams:idB": "2212",
            "Beams:eCM": "13000",
            "PartonLevel:ISR": "on",
        }

        for key, expected_val in expected.items():
            assert key in parsed, f"Key '{key}' not found"
            assert parsed[key] == expected_val


class TestPythia8CardUpdate:
    """Tests for Pythia8 card updating."""

    def test_update_pythia8_card(self, test_data_dir, tmp_path):
        """Test updating a pythia8 card with new values."""
        # Check if shower card exists (it should be pythia8 format)
        shower_card_path = test_data_dir / "shower_card.dat"
        if not shower_card_path.exists():
            pytest.skip("shower_card.dat test data not found")

        # Copy to temp directory
        test_path = tmp_path / "shower_card.dat"
        shutil.copy(shower_card_path, test_path)

        # Parse original
        original_parsed = parse_pythia8_card(test_path)

        # Update with new values
        updates = {
            "Beams:eCM": "14000",
            "Main:numberOfEvents": "500",
            # "JetMatching:qCut": "-1.0",
        }

        update_pythia8_card(test_path, updates, tmp_path)

        # Parse updated
        updated_parsed = parse_pythia8_card(test_path)

        expected = {**original_parsed, **updates}

        # Check updates were applied
        assert updated_parsed == expected

    def test_update_pythia8_card_default(self, test_data_dir, tmp_path):
        """Test updating a pythia8_card_default.dat file."""
        pythia8_card_path = test_data_dir / "pythia8_card_default.dat"
        if not pythia8_card_path.exists():
            pytest.skip("pythia8_card_default.dat test data not found")

        # Copy to temp directory
        test_path = tmp_path / "pythia8_card_default.dat"
        shutil.copy(pythia8_card_path, test_path)

        # Parse original
        original_parsed = parse_pythia8_card(test_path)

        # Update with new values - mix of existing and new parameters
        updates = {
            "Beams:idA": "2212",  # new parameter (proton)
            "Beams:idB": "2212",  # new parameter (proton)
            "Beams:eCM": "13000.0",  # new parameter (13 TeV)
            "Main:numberOfEvents": "1000",  # existing parameter (reduced number of events)
            "JetMatching:qCut": "20.0",  # existing parameter (jet matching scale)
        }

        update_pythia8_card(test_path, updates, tmp_path)

        # Parse updated
        updated_parsed = parse_pythia8_card(test_path)

        expected = {
            "Main:numberOfEvents": "1000",
            "HEPMCoutput:file": "hepmc.gz",
            "JetMatching:qCut": "20.0",
            "JetMatching:doShowerKt": "off",
            "JetMatching:nJetMax": "-1",
            "Merging:TMS": "-1.0",
            "Merging:Process": "<set_by_user>",
            "Merging:nJetMax": "-1",
            "SysCalc:fullCutVariation": "off",
            "Beams:idA": "2212",
            "Beams:idB": "2212",
            "Beams:eCM": "13000.0",
        }

        assert updated_parsed == expected
