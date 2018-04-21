#!/usr/bin/env python

from pycparser import c_ast, parse_file
from subprocess import Popen, PIPE
import sys
import os
import argparse

HERE = os.path.dirname(__file__)
DOXY_RELDIR = 'doxygen/html'


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


def parse(filename):
    ast = parse_file(filename, use_cpp=True, cpp_args=r'-Iutils/fake_libc_include')
    collector = DocItemCollector()
    collector.visit(ast)
    return collector


def generate(filename):
    d = parse(filename)

    enums = """
Constants
---------
%s

""" % '\n.. doxygenenum:: '.join([''] + d.collected['enum'])

    call_graphs = {}
    cmd = ['find', os.path.join(HERE, DOXY_RELDIR), '-regex', '.*_icgraph.dot']
    out, _ = Popen(cmd, stdout=PIPE).communicate()
    for call_graph in out.strip().split('\n'):
        with open(call_graph) as f:
            # the first line of each caller graph file is 'digraph "func"'
            func = f.readline().strip().split()[1].strip('"')
        call_graphs[func] = os.path.relpath(call_graph, HERE)

    funcs = """
Functions
---------
%s
"""
    _funcs = ''
    for func in d.collected['func']:
        _funcs += '\n\n.. doxygenfunction:: %s' % func
        if func in call_graphs:
            _funcs += '\n.. graphviz:: %s' % call_graphs[func]
    funcs = funcs % _funcs

    structs = """
Data Structures
---------------
%s
"""
    # the :align: option needs sphinx >= 1.5
    # cf. https://github.com/sphinx-doc/sphinx/commit/1d17475a
    # otherwise use the following in css:
    # img[src*="graphviz"] {
      # display: block;
      # margin: auto;
    # }
    struct_tpl = '\n' + \
            '\n.. doxygenstruct:: %s' + \
            '\n.. graphviz:: %s/struct%%s__coll__graph.dot' % DOXY_RELDIR + \
            '\n   :align: center'
    structs = structs % ''.join(struct_tpl % (s, s.replace('_', '__'))
                                for s in d.collected['struct'])

    return 'pwlib (C component)\n' + \
           '===================\n\n' + \
           'Documentation for data structrs, constants, and functions' + \
           'defined in the C component.\n' + enums + funcs + structs

def main(args):
    parser = argparse.ArgumentParser(prog='cdocs')
    parser.add_argument('source', help='The single header file to be parsed')

    namespace = parser.parse_args(args)
    sys.stdout.write(generate(namespace.source))
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
