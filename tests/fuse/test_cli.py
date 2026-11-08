"""cenfuse's command line (cli.py; future.md #65): the parsing, the error path, the control dispatch and SIGTERM."""
import logging
import os
import signal

import pytest

from cenfuse import cli
from cenfuse.util import ENGINE_LOGGERS


@pytest.fixture(autouse=True)
def _keep_logging():
    """cli.main configures the engine loggers for the process (util.setup_logging: their own handler, no propagation); put them
    back, or later tests in this worker would capture none of their records."""
    keep = {n: (logging.getLogger(n).handlers[:], logging.getLogger(n).propagate, logging.getLogger(n).level)
            for n in ENGINE_LOGGERS}
    yield
    for n, (h, prop, lv) in keep.items():
        lg = logging.getLogger(n)
        lg.handlers[:], lg.propagate = h, prop
        lg.setLevel(lv)


def test_force_takes_every_stage_and_all(capsys):
    from cenfuse.pipeline import STAGES
    from cenfuse.sat.pipeline import STAGES as SAT
    with pytest.raises(SystemExit) as e:
        cli.main(["run", "-h"])
    assert e.value.code == 0
    out = capsys.readouterr().out
    assert "--force STAGE" in out and all(s in out for s in STAGES + SAT)
    with pytest.raises(SystemExit) as e:
        cli.main(["run", "x", "--force", "nonsense"])
    assert e.value.code == 2 and "invalid choice" in capsys.readouterr().err


def test_bad_input_is_a_message_and_exit_2(tmp_path, capsys):
    rc = cli.main(["run", str(tmp_path / "no_such_folder"), "-o", str(tmp_path / "out")])
    err = capsys.readouterr().err
    assert rc == 2 and err.startswith("cenfuse: error:") and "Traceback" not in err


def test_control_is_its_own_command(capsys):
    with pytest.raises(SystemExit) as e:
        cli.main(["control", "-h"])
    assert e.value.code == 0
    out = capsys.readouterr().out
    assert "hic-repair" in out and "gc-shuffle" in out


def test_sigterm_ends_the_main_process_as_an_exception(monkeypatch):
    monkeypatch.setattr(cli, "_MAIN", os.getpid())
    with pytest.raises(SystemExit) as e:
        cli._terminate(signal.SIGTERM, None)
    assert e.value.code == 128 + signal.SIGTERM
