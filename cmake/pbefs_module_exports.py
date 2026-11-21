#!/usr/bin/env python3
###############################################################################
############################ pbefs_module_exports.py ##########################
###############################################################################
#                                                                             #
#   Writes src/PBeFS.cppm, the interface of the C++20 module pbefs, from the  #
#   AST of the headers of include/PBeFS dumped by clang: the module exports,  #
#   with a using-declaration, every name that the headers declare at the      #
#   scope of the namespace PBeFS_di_unipi_it or of one nested in it, but for  #
#   the namespaces detail:                                                    #
#                                                                             #
#       pbefs_module_exports.py [--clang <clang++>] -- -I<ADS>/include \      #
#               -DADS_WITH_MPI -I<MPI>                                        #
#                                                                             #
#   The headers are read once for each set of data structures of A* (ads, std #
#   and boost, see PBeFS/astar/Structures.h), and a name that only some of    #
#   the sets declare is exported inside the #if of their macros; the headers  #
#   that need an optional library (those of MPI, see OPTIONAL) are read only  #
#   if the flags define its macro, and their names are exported inside its    #
#   #if. To be run whenever a header gains or loses a name.                   #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Università di Pisa                              #
#                                                                             #
###############################################################################

import argparse
import collections
import json
import os
import re
import subprocess
import sys
import tempfile
import textwrap

NAMESPACE = "PBeFS_di_unipi_it"

# the header of everything that needs no optional library
ENTRY = "PBeFS/PBeFS.h"

# the headers that need an optional library, and the macro that its target
# defines [see CMakeLists.txt]; the headers they include and ENTRY does not
# need it too
OPTIONAL = {"PBeFS/astar/HDAStarMPI.h": "ADS_WITH_MPI",
            "PBeFS/astar/ParallelIDAStarMPI.h": "ADS_WITH_MPI",
            "PBeFS/bnb/MPIHubBranchAndBound.h": "ADS_WITH_MPI"}

# the sets of data structures, each with the macro that chooses it
STRUCTURES = {"ads": "PBEFS_STRUCTURES_ADS",
              "std": "PBEFS_STRUCTURES_STD",
              "boost": "PBEFS_STRUCTURES_BOOST"}

# the kinds of declaration that a using-declaration can export; the others
# (e.g., the deduction guides) come with their class
EXPORTED = {"CXXRecordDecl", "ClassTemplateDecl", "ConceptDecl", "EnumDecl",
            "FunctionDecl", "FunctionTemplateDecl", "TypedefDecl",
            "TypeAliasDecl", "TypeAliasTemplateDecl", "VarDecl",
            "VarTemplateDecl"}

INCLUDE = re.compile(r'^\s*#\s*include\s+"(PBeFS/[^"]+)"', re.M)

# the headers - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


def reach(inc, header, seen=None):
    """the headers of inc that header includes, itself among them, whatever
    the #if around the #include"""
    seen = set() if seen is None else seen
    if header in seen or not os.path.exists(os.path.join(inc, header)):
        return seen
    seen.add(header)
    with open(os.path.join(inc, header)) as f:
        for h in INCLUDE.findall(f.read()):
            reach(inc, h, seen)
    return seen


def macros(inc):
    """{header: macro} for the headers that need an optional library, and
    the headers that neither ENTRY nor these include"""
    base = reach(inc, ENTRY)
    need = {}
    for h, m in OPTIONAL.items():
        for d in reach(inc, h) - base:
            need[d] = m
    every = {os.path.relpath(os.path.join(r, f), inc)
             for r, _, fs in os.walk(os.path.join(inc, "PBeFS"))
             for f in fs if f.endswith(".h")}
    return need, sorted(every - base - set(need))

# the AST - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


def dump(clang, inc, headers, flags):
    """the declarations of the namespace PBeFS_di_unipi_it, one object of
    the JSON AST of clang for each time the headers open the namespace"""
    with tempfile.NamedTemporaryFile("w", suffix=".cpp", delete=False) as tu:
        for h in headers:
            tu.write(f'#include "{h}"\n')
    err = tempfile.TemporaryFile()
    try:
        proc = subprocess.Popen([clang, "-std=c++20", "-fsyntax-only",
                                 "-I" + inc, *flags,
                                 "-Xclang", "-ast-dump=json",
                                 "-Xclang", "-ast-dump-filter=" + NAMESPACE,
                                 tu.name],
                                stdout=subprocess.PIPE, stderr=err, text=True)
        lines = []
        for line in proc.stdout:
            lines.append(line)
            if line == "}\n":
                node = json.loads("".join(lines))
                lines = []
                if node.get("kind") == "NamespaceDecl" and \
                   node.get("name") == NAMESPACE:
                    yield node
        proc.wait()
    finally:
        os.unlink(tu.name)
    if proc.returncode:
        err.seek(0)
        errors = [ln for ln in err.read().decode().splitlines()
                  if "error" in ln]
        sys.exit("the headers do not compile:\n  " +
                 "\n  ".join(errors[:10]))


def collect(nodes, inc):
    """{(namespace, name): header} for the names that the headers of inc
    declare at namespace scope, a class under the header that defines it,
    and the names that cannot be exported"""
    names, internal, defined = {}, [], set()
    # the dump writes the file of a location only when it changes: the last
    # one is followed through every location, in the order of the dump
    state = {"file": ""}

    def follow(obj):
        if isinstance(obj, dict):
            for k, v in obj.items():
                if k == "file":
                    state["file"] = v
                elif k != "includedFrom":
                    follow(v)
        elif isinstance(obj, list):
            for v in obj:
                follow(v)

    def walk(node, ns):
        kind = node.get("kind")
        follow(node.get("loc", {}))
        path = state["file"]
        follow(node.get("range", {}))
        if kind == "NamespaceDecl":
            for n in node.get("inner", []):
                walk(n, ns + [node.get("name", "")])
            return
        follow(node.get("inner", []))
        if not os.path.abspath(path).startswith(inc + os.sep) or \
           kind not in EXPORTED or node.get("isImplicit") or \
           not node.get("name") or "detail" in ns:
            return
        name, key = node["name"], "::".join(ns)
        header = os.path.relpath(os.path.abspath(path), inc)
        # a definition of a name declared elsewhere (e.g., a member defined
        # out of its class) comes with its declaration
        if "parentDeclContextId" in node:
            return
        decl = node
        if kind == "FunctionTemplateDecl":
            decl = next((n for n in node.get("inner", [])
                         if n.get("kind", "").endswith("Decl") and
                         n.get("kind") not in ("TemplateTypeParmDecl",
                                               "NonTypeTemplateParmDecl",
                                               "TemplateTemplateParmDecl")),
                        {})
            if decl.get("kind") != "FunctionDecl":
                return
        if "" in ns or decl.get("storageClass") == "static" or \
           (kind == "VarDecl" and not node.get("inline") and
            node.get("type", {}).get("qualType", "").startswith("const")):
            internal.append(f"{key}::{name} ({header})")
            return
        record = node
        if kind == "ClassTemplateDecl":
            record = next((n for n in node.get("inner", [])
                           if n.get("kind") == "CXXRecordDecl"), {})
        if kind in ("CXXRecordDecl", "ClassTemplateDecl"):
            if record.get("completeDefinition"):
                defined.add((key, name))
                names[key, name] = header
            elif (key, name) not in defined:
                names.setdefault((key, name), None)
            return
        names.setdefault((key, name), header)
        # the enumerators of an unscoped enum live in the namespace
        if kind == "EnumDecl" and not node.get("scopedEnumTag"):
            for e in node.get("inner", []):
                if e.get("kind") == "EnumConstantDecl":
                    names.setdefault((key, e["name"]), header)

    for node in nodes:
        walk(node, [])
    # a class only declared in advance is not there
    return {k: h for k, h in names.items() if h}, internal

# the interface - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


BANNER = "/*" + "-" * 74 + "*/"


def title(text):
    text = f" {text} "
    left = (74 - len(text)) // 2
    return "/*" + "-" * left + text + "-" * (74 - left - len(text)) + "*/"


def para(text):
    return [" * " + ln for ln in textwrap.wrap(text, 75,
                                                break_on_hyphens=False)]


ABOUT = ("Interface of the C++20 module pbefs, which exports every name that "
         "the headers of include/PBeFS declare in the namespace "
         "PBeFS_di_unipi_it and in those nested in it, but for the "
         "namespaces detail; a translation unit can then write",
         "    import pbefs;",
         "in place of the #include of the headers. The names that only some "
         "sets of data structures of A* declare [see PBeFS/astar/"
         "Structures.h] are exported with the set that the module is "
         "compiled with, and those of the searches on MPI only if it is "
         "compiled with ADS_WITH_MPI, as the target PBeFS::MPI does. The "
         "macros are not exported by a module, and those that configure "
         "the headers take the value they have where the module is "
         "compiled.",
         "This file is written by cmake/pbefs_module_exports.py, which reads "
         "the names from the AST of the headers: it is not to be edited by "
         "hand.")


def condition(macros):
    return " || ".join(f"defined( {m} )" for m in macros)


def interface(found, need):
    """the lines of src/PBeFS.cppm, from {(namespace, name): {set: header}}
    and the macros of the optional headers"""
    out = [BANNER, title("File PBeFS.cppm"), BANNER, "/** @file"]
    for a in ABOUT:
        out += [" * " + a] if a.startswith(" ") else para(a)
        out.append(" *")
    out += [" * \\author Donato Meoli \\n",
            " *         Dipartimento di Informatica \\n",
            " *         Università di Pisa \\n", " *",
            " * \\copyright &copy; by Donato Meoli", " */",
            BANNER, title("GLOBAL MODULE FRAGMENT"), BANNER, "",
            "module;", "", f'#include "{ENTRY}"', ""]
    by_macro = collections.defaultdict(list)
    for h, m in need.items():
        if h in OPTIONAL:
            by_macro[m].append(h)
    for m in sorted(by_macro):
        out.append(f"#if {condition([m])}")
        out += [f' #include "{h}"' for h in sorted(by_macro[m])]
        out += ["#endif", ""]
    out += [BANNER, title("MODULE PURVIEW"), BANNER, "",
            "export module pbefs;", ""]

    # each name under the header that declares it in the first set that
    # has it, inside the #if of the sets that have it, if not all
    groups = collections.defaultdict(lambda: collections.defaultdict(list))
    for (ns, name), where in found.items():
        sets = [s for s in STRUCTURES if s in where]
        header = where[sets[0]]
        cond = []
        if len(sets) < len(STRUCTURES):
            cond.append(condition(STRUCTURES[s] for s in sets))
        if header in need:
            cond.append(condition([need[header]]))
        groups[ns][header, " && ".join(f"( {c} )" if len(cond) > 1 else c
                                         for c in cond)].append(name)
    for ns in sorted(groups):
        parts = ns.split("::")
        out += ["export namespace " + parts[0], "{"]
        for p in parts[1:]:
            out += ["namespace " + p, "{"]
        last = None
        for header, cond in sorted(groups[ns]):
            if header != last:
                out += ["", "// " + header[len("PBeFS/"):]]
                last = header
            if cond:
                out.append(f"#if {cond}")
            out += [f"using {ns}::{n};" for n in sorted(groups[ns][header,
                                                                  cond])]
            if cond:
                out.append("#endif")
        out += [""] + ["}"] * len(parts)
        out.append("")
    out += [BANNER, title("End File PBeFS.cppm"), BANNER]
    return out


def main():
    ap = argparse.ArgumentParser(
        description="Writes src/PBeFS.cppm, the interface of the module "
                    "pbefs.")
    ap.add_argument("--clang", default="clang++")
    ap.add_argument("--output", help="the file to write, src/PBeFS.cppm "
                                     "if not given")
    ap.add_argument("flags", nargs="*")
    args = ap.parse_args()

    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    inc = os.path.join(root, "include")
    need, lost = macros(inc)
    if lost:
        sys.exit("headers that neither " + ENTRY + " nor OPTIONAL include: " +
                 ", ".join(lost))
    defined = {f[2:].split("=")[0] for f in args.flags if f.startswith("-D")}
    missing = sorted(set(OPTIONAL.values()) - defined)
    if missing:
        sys.exit("the optional headers are read only with their macros, "
                 "and these are not defined: " + ", ".join(missing))
    flags = [f for f in args.flags
             if not f.startswith("-DPBEFS_STRUCTURES_")]
    headers = [ENTRY] + sorted(OPTIONAL)
    found = collections.defaultdict(dict)
    for s, m in STRUCTURES.items():
        names, internal = collect(dump(args.clang, inc, headers,
                                       flags + ["-D" + m]), inc)
        if internal:
            sys.exit("names with internal linkage, which a module cannot "
                     "export:\n  " + "\n  ".join(sorted(set(internal))))
        for k, h in names.items():
            found[k][s] = h
    out = interface(found, need)
    path = args.output or os.path.join(root, "src", "PBeFS.cppm")
    with open(path, "w") as f:
        f.write("\n".join(out) + "\n")


if __name__ == "__main__":
    main()

###############################################################################
####################### End File pbefs_module_exports.py ######################
###############################################################################
