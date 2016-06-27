#!/usr/bin/env python

# TODO documentation

from pycparser import c_ast, parse_file
import sys

def skip_name(name):
    return name[0] == '_'

class DocItemCollector(c_ast.NodeVisitor):
    def __init__(self):
        self.collected = {'enum': [], 'struct': [], 'func': []}

    # Visitor for FuncDecl nodes, collects function names
    def visit_FuncDecl(self, node):
        t = node.type
        # chase down the pointer types: e.g if return type is int** we do this
        # twice to get to the function name:
        while not isinstance(t, c_ast.TypeDecl):
            t = t.type
        if not skip_name(t.declname):
            self.collected['func'].append(t.declname)

    def visit_Typedef(self, node):
        if isinstance(node.type.type, c_ast.Enum) and not skip_name(node.name):
            self.collected['enum'].append(node.name)
        elif isinstance(node.type.type, c_ast.Struct) and not skip_name(node.name):
            self.collected['struct'].append(node.name)

def generate(filename):
    ast = parse_file(filename, use_cpp=True, cpp_args=r'-Iutils/fake_libc_include')
    d = DocItemCollector()
    d.visit(ast)
    enum = '\n.. doxygenenum:: '.join([''] + sorted(d.collected['enum']))
    func = '\n.. doxygenfunction:: '.join([''] + sorted(d.collected['func']))
    struct = '\n.. doxygenstruct:: '.join([''] + sorted(d.collected['struct']))

    enum = '\n\nConstants\n---------\n' + enum
    func = '\n\n\nFunctions\n---------\n' + func
    struct = '\n\n\nData Structures\n---------------\n' + struct

    return 'pwlib (C component)\n===================' + enum + struct + func


if __name__ == "__main__":
    sys.stdout.write(generate(sys.argv[1]))
