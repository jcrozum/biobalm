# type: ignore
import os

from biodivine_aeon import BooleanNetwork

# The purpose of this file is to detect tests with `network_file` as input and
# then supply these tests with networks from `bbm-bnet-inputs-true` up to a
# certain network size. This network size can be configured using
# `--networksize` and its default value is 20.

# We intentionally test on the `-inputs-true` models as opposed to `-inputs-identity`,
# as having fixed inputs ensures there are not too many trap spaces, fixed points, etc.


def pytest_addoption(parser):
    parser.addoption(
        "--networksize",
        action="store",
        default="20",
        help="Only check networks up to this size.",
    )


def pytest_generate_tests(metafunc):
    if "network_file" in metafunc.fixturenames:
        size = int(metafunc.config.getoption("networksize"))
        models = []
        for model in os.listdir("./models/bbm-bnet-inputs-true"):
            if not model.endswith(".bnet"):
                # Just in case there are some other files there.
                continue
            path = f"./models/bbm-bnet-inputs-true/{model}"
            bn = BooleanNetwork.from_file(path)
            if bn.num_vars() > size:
                continue
            models.append(path)
        metafunc.parametrize("network_file", models)
