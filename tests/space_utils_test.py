from biodivine_aeon import BooleanNetwork # type: ignore
from pyeda.inter import * # type: ignore

from nfvsmotifs.space_utils import *

def test_is_subspace():
    assert is_subspace({'x': '0', 'y': '1'}, {'x': '0'})
    assert not is_subspace({'x': '1', 'y': '0'}, {'x':'0', 'y':'0'})

def test_expression_percolation():
    space = {'x': '0', 'y': '1'}
    e = expr("(a & ~x) | (a & y)")    

    assert expr("a") == percolate_pyeda_expression(e, {'x': '0', 'y': '1'})
    assert expr(False) == percolate_pyeda_expression(e, {'a': '0'})

def test_space_percolation():
    bn = BooleanNetwork.from_bnet("""
        a, b
        b, c
        c, a
    """)

    assert {'a': '0', 'b': '0', 'c': '0'} == percolate_space(bn, {'a': '0'})[0]
    assert {} == percolate_space(bn, {'a': '0'})[1]
    assert is_syntactic_trap_space(bn, {'a': '0', 'b': '0', 'c': '0'})
    assert {'a': '1', 'b': '1', 'c': '1'} == percolate_space(bn, {'a': '1'})[0]
    assert {} == percolate_space(bn, {'a': '1'})[1]
    assert is_syntactic_trap_space(bn, {'a': '1', 'b': '1', 'c': '1'})

    bn = BooleanNetwork.from_bnet("""
    a, b
    b, !c
    c, a
    """)

    print(percolate_space(bn, {'a': '0', 'b': '0', 'c': '0'}))
    assert {'a': '0', 'b': '0', 'c': '0'} == percolate_space(bn, {'a': '0', 'b': '0', 'c': '0'})[0]
    # The conflict is on b---the rest is fine.
    assert {'b': '1'} == percolate_space(bn, {'a': '0', 'b': '0', 'c': '0'})[1]
    assert not is_syntactic_trap_space(bn, {'a': '0'})
    assert is_syntactic_trap_space(bn, {})

def test_network_percolation():
    bn = BooleanNetwork.from_bnet("""
        a, c & b
        b, !a
        c, c
    """)
    
    percolated_bn = percolate_network(bn, percolate_space(bn, {'c': '0'})[0])
    assert "false" == percolated_bn.get_update_function("c")
    assert "false" == percolated_bn.get_update_function("a")
    assert "true" == percolated_bn.get_update_function("b")
    percolated_bn = percolate_network(bn, percolate_space(bn, {'c': '1'})[0])
    assert "true" == percolated_bn.get_update_function("c")
    assert "b" == percolated_bn.get_update_function("a")
    assert "!a" == percolated_bn.get_update_function("b")

def test_expression_to_spaces():
    e = expr("(a & c) | (~d & (a | c)) | f")

    spaces = expression_to_space_list(e)

    assert {'f': '1'} in spaces
    assert {'a': '1', 'c': '1'} in spaces
    assert {'d': '0', 'a': '1'} in spaces
    assert {'d': '0', 'c': '1'} in spaces
    assert {'a': '1', 'c': '0'} not in spaces