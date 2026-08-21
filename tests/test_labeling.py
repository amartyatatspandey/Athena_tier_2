from secondlook.labeling import THRESHOLD_DECIDED, label_binding_score


def test_placeholder_always_returns_uncertain():
    assert THRESHOLD_DECIDED is False
    assert label_binding_score(-2.056) == "uncertain"
    assert label_binding_score(0.0) == "uncertain"
    assert label_binding_score(1.5) == "uncertain"
    assert label_binding_score(None) == "uncertain"
