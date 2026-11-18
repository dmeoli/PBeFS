begin_version
3
end_version
begin_metric
0
end_metric
9
begin_variable
var0
-1
6
Atom clear(a)
Atom holding(a)
Atom on(a, a)
Atom on(b, a)
Atom on(c, a)
Atom on(d, a)
end_variable
begin_variable
var1
-1
6
Atom clear(b)
Atom holding(b)
Atom on(a, b)
Atom on(b, b)
Atom on(c, b)
Atom on(d, b)
end_variable
begin_variable
var2
-1
6
Atom clear(c)
Atom holding(c)
Atom on(a, c)
Atom on(b, c)
Atom on(c, c)
Atom on(d, c)
end_variable
begin_variable
var3
-1
6
Atom clear(d)
Atom holding(d)
Atom on(a, d)
Atom on(b, d)
Atom on(c, d)
Atom on(d, d)
end_variable
begin_variable
var4
-1
2
Atom handempty()
NegatedAtom handempty()
end_variable
begin_variable
var5
-1
2
Atom ontable(a)
NegatedAtom ontable(a)
end_variable
begin_variable
var6
-1
2
Atom ontable(b)
NegatedAtom ontable(b)
end_variable
begin_variable
var7
-1
2
Atom ontable(c)
NegatedAtom ontable(c)
end_variable
begin_variable
var8
-1
2
Atom ontable(d)
NegatedAtom ontable(d)
end_variable
4
begin_mutex_group
6
0 0
0 1
0 2
0 3
0 4
0 5
end_mutex_group
begin_mutex_group
6
1 0
1 1
1 2
1 3
1 4
1 5
end_mutex_group
begin_mutex_group
6
2 0
2 1
2 2
2 3
2 4
2 5
end_mutex_group
begin_mutex_group
6
3 0
3 1
3 2
3 3
3 4
3 5
end_mutex_group
begin_state
0
0
0
0
0
0
0
0
0
end_state
begin_goal
3
0 3
1 4
2 5
end_goal
32
begin_operator
pick-up a
0
3
0 0 0 1
0 4 0 1
0 5 0 1
1
end_operator
begin_operator
pick-up b
0
3
0 1 0 1
0 4 0 1
0 6 0 1
1
end_operator
begin_operator
pick-up c
0
3
0 2 0 1
0 4 0 1
0 7 0 1
1
end_operator
begin_operator
pick-up d
0
3
0 3 0 1
0 4 0 1
0 8 0 1
1
end_operator
begin_operator
put-down a
0
3
0 0 1 0
0 4 -1 0
0 5 -1 0
1
end_operator
begin_operator
put-down b
0
3
0 1 1 0
0 4 -1 0
0 6 -1 0
1
end_operator
begin_operator
put-down c
0
3
0 2 1 0
0 4 -1 0
0 7 -1 0
1
end_operator
begin_operator
put-down d
0
3
0 3 1 0
0 4 -1 0
0 8 -1 0
1
end_operator
begin_operator
stack a b
0
3
0 0 1 0
0 1 0 2
0 4 -1 0
1
end_operator
begin_operator
stack a c
0
3
0 0 1 0
0 2 0 2
0 4 -1 0
1
end_operator
begin_operator
stack a d
0
3
0 0 1 0
0 3 0 2
0 4 -1 0
1
end_operator
begin_operator
stack b a
0
3
0 0 0 3
0 1 1 0
0 4 -1 0
1
end_operator
begin_operator
stack b c
0
3
0 1 1 0
0 2 0 3
0 4 -1 0
1
end_operator
begin_operator
stack b d
0
3
0 1 1 0
0 3 0 3
0 4 -1 0
1
end_operator
begin_operator
stack c a
0
3
0 0 0 4
0 2 1 0
0 4 -1 0
1
end_operator
begin_operator
stack c b
0
3
0 1 0 4
0 2 1 0
0 4 -1 0
1
end_operator
begin_operator
stack c d
0
3
0 2 1 0
0 3 0 4
0 4 -1 0
1
end_operator
begin_operator
stack d a
0
3
0 0 0 5
0 3 1 0
0 4 -1 0
1
end_operator
begin_operator
stack d b
0
3
0 1 0 5
0 3 1 0
0 4 -1 0
1
end_operator
begin_operator
stack d c
0
3
0 2 0 5
0 3 1 0
0 4 -1 0
1
end_operator
begin_operator
unstack a b
0
3
0 0 0 1
0 1 2 0
0 4 0 1
1
end_operator
begin_operator
unstack a c
0
3
0 0 0 1
0 2 2 0
0 4 0 1
1
end_operator
begin_operator
unstack a d
0
3
0 0 0 1
0 3 2 0
0 4 0 1
1
end_operator
begin_operator
unstack b a
0
3
0 0 3 0
0 1 0 1
0 4 0 1
1
end_operator
begin_operator
unstack b c
0
3
0 1 0 1
0 2 3 0
0 4 0 1
1
end_operator
begin_operator
unstack b d
0
3
0 1 0 1
0 3 3 0
0 4 0 1
1
end_operator
begin_operator
unstack c a
0
3
0 0 4 0
0 2 0 1
0 4 0 1
1
end_operator
begin_operator
unstack c b
0
3
0 1 4 0
0 2 0 1
0 4 0 1
1
end_operator
begin_operator
unstack c d
0
3
0 2 0 1
0 3 4 0
0 4 0 1
1
end_operator
begin_operator
unstack d a
0
3
0 0 5 0
0 3 0 1
0 4 0 1
1
end_operator
begin_operator
unstack d b
0
3
0 1 5 0
0 3 0 1
0 4 0 1
1
end_operator
begin_operator
unstack d c
0
3
0 2 5 0
0 3 0 1
0 4 0 1
1
end_operator
0
