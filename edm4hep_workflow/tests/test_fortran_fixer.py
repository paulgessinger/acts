"""Tests for Fortran $ syntax fixer."""

import pytest
from pathlib import Path
import shutil
from src.colliderml.cli.madgraph import fix_fortran_dollar_syntax


class TestFortranFixer:
    """Tests for Fortran $ syntax fixer."""

    def test_fix_fortran_dollar_syntax(self, tmp_path):
        """Test fixing Fortran $ format descriptors."""
        # Create test file with $ syntax
        test_content = """      subroutine test()
            if (mapconfig(i) .lt. 10) then
               write(26,'(x,i1,a2$)') mapconfig(i),postfix
            elseif (mapconfig(i) .lt. 100) then
               write(26,'(x,i2,a2$)') mapconfig(i),postfix
            endif
      write(lun,'(a$)') 'for i in $channel '
      end subroutine"""

        test_file = tmp_path / "test_fortran.f"
        test_file.write_text(test_content)

        # Apply fix
        fix_fortran_dollar_syntax(test_file, tmp_path)

        # Read fixed content
        fixed_content = test_file.read_text()

        # Check that $ syntax was replaced
        assert "$)'" not in fixed_content
        assert "advance='no'" in fixed_content

        # Check specific replacements
        assert (
            "write(26,'(x,i1,a2)',advance='no') mapconfig(i),postfix" in fixed_content
        )
        assert (
            "write(26,'(x,i2,a2)',advance='no') mapconfig(i),postfix" in fixed_content
        )
        assert "write(lun,'(a)',advance='no') 'for i in $channel '" in fixed_content

    def test_fix_fortran_no_changes_needed(self, tmp_path):
        """Test that files without $ syntax are not modified."""
        # Create test file without $ syntax
        test_content = """      subroutine test()
      write(26,'(a)') 'normal fortran'
      write(lun,'(i0)',advance='no') 42
      end subroutine"""

        test_file = tmp_path / "test_normal.f"
        test_file.write_text(test_content)
        original_content = test_content

        # Apply fix (should do nothing)
        fix_fortran_dollar_syntax(test_file, tmp_path)

        # Content should be unchanged
        fixed_content = test_file.read_text()
        assert fixed_content == original_content

    def test_fix_fortran_file_not_found(self, tmp_path):
        """Test handling of non-existent files."""
        non_existent = tmp_path / "does_not_exist.f"

        # Should not raise exception
        fix_fortran_dollar_syntax(non_existent, tmp_path)

    def test_fix_multiple_dollar_formats(self, tmp_path):
        """Test fixing multiple different $ format variations."""
        test_content = """      program test
      write(10,'(i5$)') num1
      write(20,'(f10.3$)') real_val
      write(30,'(a10$)') text_val
      write(40,'(e12.4$)') sci_val
      end program"""

        test_file = tmp_path / "test_multiple.f"
        test_file.write_text(test_content)

        # Apply fix
        fix_fortran_dollar_syntax(test_file, tmp_path)

        # Read fixed content
        fixed_content = test_file.read_text()

        # Check all variations were fixed
        expected_replacements = [
            "write(10,'(i5)',advance='no') num1",
            "write(20,'(f10.3)',advance='no') real_val", 
            "write(30,'(a10)',advance='no') text_val",
            "write(40,'(e12.4)',advance='no') sci_val"
        ]

        for expected in expected_replacements:
            assert expected in fixed_content

        # Ensure no $ syntax remains
        assert "$)'" not in fixed_content

    def test_fix_complex_format_strings(self, tmp_path):
        """Test fixing complex format strings with multiple specifiers."""
        test_content = """      subroutine complex_output()
      write(unit,'(i3,1x,f8.2,1x,a,$)') id, value, label, extra
      write(26,'(2i4,3f10.3$)') n1, n2, x, y, z
      end subroutine"""

        test_file = tmp_path / "test_complex.f"
        test_file.write_text(test_content)

        # Apply fix
        fix_fortran_dollar_syntax(test_file, tmp_path)

        # Read fixed content
        fixed_content = test_file.read_text()

        # Check complex formats were handled correctly
        assert "write(unit,'(i3,1x,f8.2,1x,a,)',advance='no') id, value, label, extra" in fixed_content
        assert "write(26,'(2i4,3f10.3)',advance='no') n1, n2, x, y, z" in fixed_content
        assert "$)'" not in fixed_content