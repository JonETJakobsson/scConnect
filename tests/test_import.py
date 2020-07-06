import scConnect as cn
import pytest

def test_import():
    assert isinstance(type(cn), "module")
    