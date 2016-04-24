import pytest

# working example of py.test usage with a fixture

class Foo():
    def __init__(self, i):
        self.i = i

@pytest.fixture(params=[Foo(i) for i in range(3,7)])
def alnprob(request):
    return request.param

def test_foo(alnprob):
    assert alnprob.i == 5
