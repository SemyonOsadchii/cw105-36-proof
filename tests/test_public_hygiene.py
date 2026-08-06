from __future__ import annotations

from tools.check_public_hygiene import ROOT, audit


def test_tracked_tree_has_no_private_agent_artifacts() -> None:
    assert audit(ROOT) == []
