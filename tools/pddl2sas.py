#!/usr/bin/env python3
###############################################################################
################################# pddl2sas.py #################################
###############################################################################
#                                                                             #
#   Translates a STRIPS task in PDDL (with :typing, without negative          #
#   preconditions, conditional effects or quantifiers) into the output.sas    #
#   format of the Fast Downward translator, version 3:                        #
#                                                                             #
#       pddl2sas.py domain.pddl problem.pddl output.sas                       #
#                                                                             #
#   The actions are grounded by relaxed reachability, and the fluent atoms    #
#   are grouped into finite-domain variables by the invariants "at most one   #
#   atom of the group holds" that are proved on the ground actions; each      #
#   atom left out of the groups becomes a binary variable.                    #
#                                                                             #
#                                Donato Meoli                                 #
#                         Dipartimento di Informatica                         #
#                             Università di Pisa                              #
#                                                                             #
###############################################################################

import itertools
import sys

# the s-expressions - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


def parse_sexp(text):
    """the nested lists of the lowercase tokens of text, comments removed"""
    lines = [ln.split(";", 1)[0] for ln in text.lower().splitlines()]
    tokens = " ".join(lines).replace("(", " ( ").replace(")", " ) ").split()
    stack = [[]]
    for t in tokens:
        if t == "(":
            stack.append([])
        elif t == ")":
            done = stack.pop()
            stack[-1].append(done)
        else:
            stack[-1].append(t)
    return stack[0][0]


def typed_list(items):
    """the pairs (name, type) of a typed list, untyped names of type object"""
    out, names, i = [], [], 0
    while i < len(items):
        if items[i] == "-":
            typ = items[i + 1]
            out += [(n, typ if isinstance(typ, str) else "object")
                    for n in names]
            names, i = [], i + 2
        else:
            names.append(items[i])
            i += 1
    return out + [(n, "object") for n in names]


def conjunction(e):
    """the literals of a conjunction, each as (positive, atom)"""
    if not e:
        return []
    if e[0] == "and":
        return [x for c in e[1:] for x in conjunction(c)]
    if e[0] == "not":
        return [(False, tuple(e[1]))]
    if e[0] in ("or", "imply", "forall", "exists", "when"):
        raise SystemExit("pddl2sas: only STRIPS tasks, not " + e[0])
    return [(True, tuple(e))]


# the task - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


def read_task(dom_text, prob_text):
    dom, prob = parse_sexp(dom_text), parse_sexp(prob_text)
    parents, consts, actions = {}, [], []
    for sec in dom[2:]:
        if sec[0] == ":types":
            parents.update(typed_list(sec[1:]))
        elif sec[0] == ":constants":
            consts += typed_list(sec[1:])
        elif sec[0] == ":action":
            a = {"name": sec[1], "params": [], "pre": [], "eff": []}
            for k in range(2, len(sec), 2):
                if sec[k] == ":parameters":
                    a["params"] = typed_list(sec[k + 1])
                elif sec[k] == ":precondition":
                    a["pre"] = conjunction(sec[k + 1])
                elif sec[k] == ":effect":
                    a["eff"] = conjunction(sec[k + 1])
            if any(not p for p, _ in a["pre"]):
                raise SystemExit("pddl2sas: negative preconditions")
            actions.append(a)
    objs, init, goal = list(consts), [], []
    for sec in prob[2:]:
        if sec[0] == ":objects":
            objs += typed_list(sec[1:])
        elif sec[0] == ":init":
            init = [tuple(x) for x in sec[1:]]
        elif sec[0] == ":goal":
            goal = conjunction(sec[1])
    if any(not p for p, _ in goal):
        raise SystemExit("pddl2sas: negative goals")

    def ancestors(t):
        while t in parents and parents[t] != t:
            yield t
            t = parents[t]
        yield t

    by_type = {}
    for o, t in objs:
        for a in ancestors(t):
            by_type.setdefault(a, []).append(o)
        by_type.setdefault("object", []).append(o)
    for t in by_type:
        by_type[t] = sorted(set(by_type[t]))
    return actions, by_type, set(init), [a for _, a in goal]


# grounding by relaxed reachability - - - - - - - - - - - - - - - - - - - - -


def ground(actions, by_type, init):
    """the ground actions reachable in the relaxed task and the atoms"""
    reach = set(init)
    found = {}
    changed = True
    while changed:
        changed = False
        by_pred = {}
        for at in reach:
            by_pred.setdefault(at[0], []).append(at)
        for a in actions:
            for b in bindings(a, by_type, by_pred):
                key = (a["name"],) + tuple(b[p] for p, _ in a["params"])
                if key in found:
                    continue
                inst = lambda at: (at[0],) + tuple(b.get(x, x)
                                                   for x in at[1:])
                pre = sorted({inst(at) for _, at in a["pre"]
                              if at[0] != "="})
                add = sorted({inst(at) for p, at in a["eff"] if p})
                dele = sorted({inst(at) for p, at in a["eff"] if not p}
                              - set(add))
                found[key] = (pre, add, dele)
                for at in add:
                    if at not in reach:
                        reach.add(at)
                        changed = True
    return found, reach


def bindings(a, by_type, by_pred):
    """the bindings of the parameters of a satisfying its preconditions in
    the atoms by_pred, by a join over the preconditions"""
    types = dict(a["params"])
    pre = [at for _, at in a["pre"] if at[0] != "="]
    eqs = [at for _, at in a["pre"] if at[0] == "="]
    pre.sort(key=lambda at: len(by_pred.get(at[0], [])))
    out = []

    def rec(k, b):
        if k == len(pre):
            free = [p for p in types if p not in b]
            for vals in itertools.product(*[by_type.get(types[p], [])
                                            for p in free]):
                bb = dict(b)
                bb.update(zip(free, vals))
                if all(bb.get(x, x) == bb.get(y, y) for _, x, y in eqs):
                    out.append(bb)
            return
        at = pre[k]
        for cand in by_pred.get(at[0], []):
            if len(cand) != len(at):
                continue
            bb, ok = dict(b), True
            for x, v in zip(at[1:], cand[1:]):
                if x.startswith("?"):
                    if bb.setdefault(x, v) != v or \
                       v not in by_type.get(types[x], []):
                        ok = False
                        break
                elif x != v:
                    ok = False
                    break
            if ok:
                rec(k + 1, bb)

    rec(0, {})
    return out


# the invariants - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


def candidate_parts(fluents, arity):
    """the parts (predicate, positions of the parameters of the group)"""
    parts = []
    for p in sorted(fluents):
        for k in range(0, min(arity[p], 2) + 1):
            for pos in itertools.permutations(range(arity[p]), k):
                parts.append((p, pos))
    return parts


def group_of(atom, invariant):
    """the key of the group of the invariant that atom belongs to, or None"""
    for p, pos in invariant:
        if atom[0] == p:
            return tuple(atom[1 + i] for i in pos)
    return None


def conflicting(pre, invariant):
    """true if pre holds two atoms of a group of the invariant"""
    pg = {}
    for at in pre:
        g = group_of(at, invariant)
        if g is not None and pg.setdefault(g, at) != at:
            return True
    return False


def proved(invariant, ops, init):
    """true if no state reachable from init has two atoms of a group, and
    each action deletes an atom of a group only if it is known to hold"""
    seen = set()
    for at in init:
        g = group_of(at, invariant)
        if g is not None:
            if g in seen:
                return False
            seen.add(g)
    for pre, add, dele in ops:
        if conflicting(pre, invariant):
            continue                           # never applicable
        pg = {}
        for at in pre:
            g = group_of(at, invariant)
            if g is not None:
                pg[g] = at
        ag = set()
        for at in add:
            g = group_of(at, invariant)
            if g is None:
                continue
            if g in ag or g not in pg:
                return False
            ag.add(g)
            if pg[g] != at and pg[g] not in dele:
                return False
        for at in dele:
            g = group_of(at, invariant)
            if g is not None and g not in pg:
                return False
    return True


def invariants(fluents, arity, ops, init):
    parts = candidate_parts(fluents, arity)
    found = []
    for n in (1, 2, 3):
        for combo in itertools.combinations(parts, n):
            k = len(combo[0][1])
            if any(len(pos) != k for _, pos in combo):
                continue
            if len({p for p, _ in combo}) != n:
                continue
            if n == 1 and k == arity[combo[0][0]]:
                continue                       # a single atom per group
            if proved(list(combo), ops, init):
                found.append(list(combo))
    return found


# the output - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


def atom_name(at):
    return "%s(%s)" % (at[0], ", ".join(at[1:]))


def translate(dom_text, prob_text):
    actions, by_type, init, goal = read_task(dom_text, prob_text)
    found, reach = ground(actions, by_type, init)
    fluents = {p for a in actions for _, at in a["eff"] for p in [at[0]]}
    arity = {}
    for at in reach | set(goal):
        arity[at[0]] = len(at) - 1
    atoms = sorted(at for at in reach if at[0] in fluents)
    ops = list(found.values())
    inv = invariants(fluents & set(arity), arity, ops, init)
    found = {k: o for k, o in found.items()
             if not any(conflicting(o[0], i) for i in inv)}

    # the groups, the largest first, each atom in one variable only
    groups = {}
    for i, invariant in enumerate(inv):
        for at in atoms:
            g = group_of(at, invariant)
            if g is not None:
                groups.setdefault((i, g), []).append(at)
    covered, variables = set(), []
    for key in sorted(groups, key=lambda x: (-len(groups[x]), x)):
        rest = [at for at in groups[key] if at not in covered]
        if len(rest) >= 2:
            variables.append(rest)
            covered.update(rest)
    for at in atoms:
        if at not in covered:
            variables.append([at])
    var_of = {}
    for v, group in enumerate(variables):
        for i, at in enumerate(group):
            var_of[at] = (v, i)
    # the operators, NONE standing for the value "none of those"
    NONE = -2
    sas_ops = []
    for key in sorted(found):
        pre, add, dele = found[key]
        if any(at not in var_of and at not in init for at in pre):
            continue                           # a static atom that is false
        pv = {}
        ok = True
        for at in pre:
            if at in var_of:
                v, i = var_of[at]
                if pv.setdefault(v, i) != i:
                    ok = False
        if not ok:
            continue
        eff = {}
        for at in dele:
            v, i = var_of[at]
            if pv.get(v, i) == i:
                eff[v] = NONE
        for at in add:
            v, i = var_of[at]
            eff[v] = i
        prevail = sorted((v, i) for v, i in pv.items() if v not in eff)
        prepost = sorted((v, pv.get(v, -1), i) for v, i in eff.items()
                         if pv.get(v, -1) != i)
        if not prepost:
            continue
        lines = ["begin_operator", " ".join(key), str(len(prevail))]
        lines += ["%d %d" % x for x in prevail]
        lines.append(str(len(prepost)))
        lines += [prepost]
        lines += ["1", "end_operator"]
        sas_ops.append(lines)
    # the value "none of those" only for a variable that may lose its atoms
    need = [not any(at in init for at in g) for g in variables]
    for lines in sas_ops:
        for v, _, post in lines[-3]:
            if post == NONE:
                need[v] = True
    none = [len(g) for g in variables]
    for lines in sas_ops:
        lines[-3] = "\n".join("0 %d %d %d" % (v, p, none[v] if q == NONE
                                                else q)
                              for v, p, q in lines[-3])

    out = ["begin_version", "3", "end_version", "begin_metric", "0",
           "end_metric", str(len(variables))]
    for v, group in enumerate(variables):
        out += ["begin_variable", "var%d" % v, "-1",
                str(len(group) + int(need[v] or len(group) == 1))]
        out += ["Atom " + atom_name(at) for at in group]
        if len(group) == 1:
            out.append("NegatedAtom " + atom_name(group[0]))
        elif need[v]:
            out.append("<none of those>")
        out.append("end_variable")
    multi = [g for g in variables if len(g) > 1]
    out.append(str(len(multi)))
    for g in multi:
        out += ["begin_mutex_group", str(len(g))]
        out += ["%d %d" % var_of[at] for at in g]
        out.append("end_mutex_group")
    out.append("begin_state")
    for v, group in enumerate(variables):
        val = [i for i, at in enumerate(group) if at in init]
        out.append(str(val[0] if val else none[v]))
    out.append("end_state")
    gl = sorted({var_of[at] for at in goal if at in var_of})
    if len(gl) != len({at for at in goal}):
        raise SystemExit("pddl2sas: a goal is unreachable or static")
    out += ["begin_goal", str(len(gl))] + ["%d %d" % x for x in gl]
    out.append("end_goal")

    out.append(str(len(sas_ops)))
    for lines in sas_ops:
        out += lines
    out.append("0")
    return "\n".join(out) + "\n"


if __name__ == "__main__":
    if len(sys.argv) != 4:
        raise SystemExit("usage: pddl2sas.py domain.pddl problem.pddl "
                         "output.sas")
    with open(sys.argv[1]) as d, open(sys.argv[2]) as p:
        text = translate(d.read(), p.read())
    with open(sys.argv[3], "w") as o:
        o.write(text)
