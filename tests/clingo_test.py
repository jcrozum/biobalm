import subprocess


def test_clingo():
    """
    This just verifies that we have clingo installed in a way that
    Trappist needs.
    """
    result = subprocess.run(
        [
            "clingo",
            "--version",
        ],
        capture_output=False,
        text=True,
    )

    assert result.returncode == 0
