"""Suite-wide: CENFUSE_REQUIRE_TOOLS=1 (CI, releases) turns every skip for a missing external tool into a failure, so a
machine without the bioconda tools cannot report green while the satellite engine's end-to-end tests never ran (future.md #65)."""
import os

# as cenfuse/__init__.py, before anything imports numpy: on a host with fragmented memory, transparent huge pages make large
# sorts and allocations many times slower (np.unique over a synthetic genome's K-mers: minutes instead of seconds)
os.environ.setdefault("NUMPY_MADVISE_HUGEPAGE", "0")

import pytest  # noqa: E402

import re

REQUIRE_TOOLS = os.environ.get("CENFUSE_REQUIRE_TOOLS") == "1"
# the suite's skip reasons for an absent external tool ("blast missing", "FasTAN tools missing", "minimap2 not on PATH", ...)
TOOL_SKIP = re.compile(r"\b(tools? missing|missing|not on PATH|not installed|not found)\b", re.IGNORECASE)


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_makereport(item, call):
    outcome = yield
    rep = outcome.get_result()
    if REQUIRE_TOOLS and rep.skipped and TOOL_SKIP.search(str(rep.longrepr)):
        rep.outcome = "failed"
        rep.longrepr = f"CENFUSE_REQUIRE_TOOLS=1: this test was skipped for a missing tool: {rep.longrepr}"
