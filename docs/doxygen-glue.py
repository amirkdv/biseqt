#!/usr/bin/env python

from pycparser import c_ast, parse_file
import sys
import argparse


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
    ret = """
pwlib (C component)
===================

Data Structures
---------------
%s

Functions
---------
%s

Constants
---------
%s
"""
    enum = '\n'.join('.. doxygenenum:: %s' % e
                     for e in d.collected['enum'])
    func = '\n'.join('.. doxygenfunction:: %s' % f
                     for f in d.collected['func'])

    # TODO with the next relase of sphinx (1.5) we can center align graphviz
    # cf. https://github.com/sphinx-doc/sphinx/commit/1d17475a
    # currently this is hacked to work in CSS, cf. docs/_static/theme_hacks.css
    struct_tpl = """
.. doxygenstruct:: %s
.. graphviz:: doxygen/html/struct%s__coll__graph.dot

---------------
"""
    struct = ''
    for s in d.collected['struct']:
        struct += (struct_tpl % (s, s.replace('_', '__')))

    return ret % (struct, func, enum)

def main(args):
    parser = argparse.ArgumentParser(prog='cdocs')
    parser.add_argument('source', help='The single header file to be parsed')

    namespace = parser.parse_args(args)
    sys.stdout.write(generate(namespace.source))
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
