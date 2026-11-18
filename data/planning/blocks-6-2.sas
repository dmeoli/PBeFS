begin_version
3
end_version
begin_metric
0
end_metric
13
begin_variable
var0
-1
8
Atom clear(a)
Atom holding(a)
Atom on(a, a)
Atom on(b, a)
Atom on(c, a)
Atom on(d, a)
Atom on(e, a)
Atom on(f, a)
end_variable
begin_variable
var1
-1
8
Atom clear(b)
Atom holding(b)
Atom on(a, b)
Atom on(b, b)
Atom on(c, b)
Atom on(d, b)
Atom on(e, b)
Atom on(f, b)
end_variable
begin_variable
var2
-1
8
Atom clear(c)
Atom holding(c)
Atom on(a, c)
Atom on(b, c)
Atom on(c, c)
Atom on(d, c)
Atom on(e, c)
Atom on(f, c)
end_variable
begin_variable
var3
-1
8
Atom clear(d)
Atom holding(d)
Atom on(a, d)
Atom on(b, d)
Atom on(c, d)
Atom on(d, d)
Atom on(e, d)
Atom on(f, d)
end_variable
begin_variable
var4
-1
8
Atom clear(e)
Atom holding(e)
Atom on(a, e)
Atom on(b, e)
Atom on(c, e)
Atom on(d, e)
Atom on(e, e)
Atom on(f, e)
end_variable
begin_variable
var5
-1
8
Atom clear(f)
Atom holding(f)
Atom on(a, f)
Atom on(b, f)
Atom on(c, f)
Atom on(d, f)
Atom on(e, f)
Atom on(f, f)
end_variable
begin_variable
var6
-1
2
Atom handempty()
NegatedAtom handempty()
end_variable
begin_variable
var7
-1
2
Atom ontable(a)
NegatedAtom ontable(a)
end_variable
begin_variable
var8
-1
2
Atom ontable(b)
NegatedAtom ontable(b)
end_variable
begin_variable
var9
-1
2
Atom ontable(c)
NegatedAtom ontable(c)
end_variable
begin_variable
var10
-1
2
Atom ontable(d)
NegatedAtom ontable(d)
end_variable
begin_variable
var11
-1
2
Atom ontable(e)
NegatedAtom ontable(e)
end_variable
begin_variable
var12
-1
2
Atom ontable(f)
NegatedAtom ontable(f)
end_variable
6
begin_mutex_group
8
0 0
0 1
0 2
0 3
0 4
0 5
0 6
0 7
end_mutex_group
begin_mutex_group
8
1 0
1 1
1 2
1 3
1 4
1 5
1 6
1 7
end_mutex_group
begin_mutex_group
8
2 0
2 1
2 2
2 3
2 4
2 5
2 6
2 7
end_mutex_group
begin_mutex_group
8
3 0
3 1
3 2
3 3
3 4
3 5
3 6
3 7
end_mutex_group
begin_mutex_group
8
4 0
4 1
4 2
4 3
4 4
4 5
4 6
4 7
end_mutex_group
begin_mutex_group
8
5 0
5 1
5 2
5 3
5 4
5 5
5 6
5 7
end_mutex_group
begin_state
0
5
6
2
7
3
0
1
1
0
1
1
1
end_state
begin_goal
5
0 7
1 2
2 3
3 4
5 6
end_goal
72
begin_operator
pick-up a
0
3
0 0 0 1
0 6 0 1
0 7 0 1
1
end_operator
begin_operator
pick-up b
0
3
0 1 0 1
0 6 0 1
0 8 0 1
1
end_operator
begin_operator
pick-up c
0
3
0 2 0 1
0 6 0 1
0 9 0 1
1
end_operator
begin_operator
pick-up d
0
3
0 3 0 1
0 6 0 1
0 10 0 1
1
end_operator
begin_operator
pick-up e
0
3
0 4 0 1
0 6 0 1
0 11 0 1
1
end_operator
begin_operator
pick-up f
0
3
0 5 0 1
0 6 0 1
0 12 0 1
1
end_operator
begin_operator
put-down a
0
3
0 0 1 0
0 6 -1 0
0 7 -1 0
1
end_operator
begin_operator
put-down b
0
3
0 1 1 0
0 6 -1 0
0 8 -1 0
1
end_operator
begin_operator
put-down c
0
3
0 2 1 0
0 6 -1 0
0 9 -1 0
1
end_operator
begin_operator
put-down d
0
3
0 3 1 0
0 6 -1 0
0 10 -1 0
1
end_operator
begin_operator
put-down e
0
3
0 4 1 0
0 6 -1 0
0 11 -1 0
1
end_operator
begin_operator
put-down f
0
3
0 5 1 0
0 6 -1 0
0 12 -1 0
1
end_operator
begin_operator
stack a b
0
3
0 0 1 0
0 1 0 2
0 6 -1 0
1
end_operator
begin_operator
stack a c
0
3
0 0 1 0
0 2 0 2
0 6 -1 0
1
end_operator
begin_operator
stack a d
0
3
0 0 1 0
0 3 0 2
0 6 -1 0
1
end_operator
begin_operator
stack a e
0
3
0 0 1 0
0 4 0 2
0 6 -1 0
1
end_operator
begin_operator
stack a f
0
3
0 0 1 0
0 5 0 2
0 6 -1 0
1
end_operator
begin_operator
stack b a
0
3
0 0 0 3
0 1 1 0
0 6 -1 0
1
end_operator
begin_operator
stack b c
0
3
0 1 1 0
0 2 0 3
0 6 -1 0
1
end_operator
begin_operator
stack b d
0
3
0 1 1 0
0 3 0 3
0 6 -1 0
1
end_operator
begin_operator
stack b e
0
3
0 1 1 0
0 4 0 3
0 6 -1 0
1
end_operator
begin_operator
stack b f
0
3
0 1 1 0
0 5 0 3
0 6 -1 0
1
end_operator
begin_operator
stack c a
0
3
0 0 0 4
0 2 1 0
0 6 -1 0
1
end_operator
begin_operator
stack c b
0
3
0 1 0 4
0 2 1 0
0 6 -1 0
1
end_operator
begin_operator
stack c d
0
3
0 2 1 0
0 3 0 4
0 6 -1 0
1
end_operator
begin_operator
stack c e
0
3
0 2 1 0
0 4 0 4
0 6 -1 0
1
end_operator
begin_operator
stack c f
0
3
0 2 1 0
0 5 0 4
0 6 -1 0
1
end_operator
begin_operator
stack d a
0
3
0 0 0 5
0 3 1 0
0 6 -1 0
1
end_operator
begin_operator
stack d b
0
3
0 1 0 5
0 3 1 0
0 6 -1 0
1
end_operator
begin_operator
stack d c
0
3
0 2 0 5
0 3 1 0
0 6 -1 0
1
end_operator
begin_operator
stack d e
0
3
0 3 1 0
0 4 0 5
0 6 -1 0
1
end_operator
begin_operator
stack d f
0
3
0 3 1 0
0 5 0 5
0 6 -1 0
1
end_operator
begin_operator
stack e a
0
3
0 0 0 6
0 4 1 0
0 6 -1 0
1
end_operator
begin_operator
stack e b
0
3
0 1 0 6
0 4 1 0
0 6 -1 0
1
end_operator
begin_operator
stack e c
0
3
0 2 0 6
0 4 1 0
0 6 -1 0
1
end_operator
begin_operator
stack e d
0
3
0 3 0 6
0 4 1 0
0 6 -1 0
1
end_operator
begin_operator
stack e f
0
3
0 4 1 0
0 5 0 6
0 6 -1 0
1
end_operator
begin_operator
stack f a
0
3
0 0 0 7
0 5 1 0
0 6 -1 0
1
end_operator
begin_operator
stack f b
0
3
0 1 0 7
0 5 1 0
0 6 -1 0
1
end_operator
begin_operator
stack f c
0
3
0 2 0 7
0 5 1 0
0 6 -1 0
1
end_operator
begin_operator
stack f d
0
3
0 3 0 7
0 5 1 0
0 6 -1 0
1
end_operator
begin_operator
stack f e
0
3
0 4 0 7
0 5 1 0
0 6 -1 0
1
end_operator
begin_operator
unstack a b
0
3
0 0 0 1
0 1 2 0
0 6 0 1
1
end_operator
begin_operator
unstack a c
0
3
0 0 0 1
0 2 2 0
0 6 0 1
1
end_operator
begin_operator
unstack a d
0
3
0 0 0 1
0 3 2 0
0 6 0 1
1
end_operator
begin_operator
unstack a e
0
3
0 0 0 1
0 4 2 0
0 6 0 1
1
end_operator
begin_operator
unstack a f
0
3
0 0 0 1
0 5 2 0
0 6 0 1
1
end_operator
begin_operator
unstack b a
0
3
0 0 3 0
0 1 0 1
0 6 0 1
1
end_operator
begin_operator
unstack b c
0
3
0 1 0 1
0 2 3 0
0 6 0 1
1
end_operator
begin_operator
unstack b d
0
3
0 1 0 1
0 3 3 0
0 6 0 1
1
end_operator
begin_operator
unstack b e
0
3
0 1 0 1
0 4 3 0
0 6 0 1
1
end_operator
begin_operator
unstack b f
0
3
0 1 0 1
0 5 3 0
0 6 0 1
1
end_operator
begin_operator
unstack c a
0
3
0 0 4 0
0 2 0 1
0 6 0 1
1
end_operator
begin_operator
unstack c b
0
3
0 1 4 0
0 2 0 1
0 6 0 1
1
end_operator
begin_operator
unstack c d
0
3
0 2 0 1
0 3 4 0
0 6 0 1
1
end_operator
begin_operator
unstack c e
0
3
0 2 0 1
0 4 4 0
0 6 0 1
1
end_operator
begin_operator
unstack c f
0
3
0 2 0 1
0 5 4 0
0 6 0 1
1
end_operator
begin_operator
unstack d a
0
3
0 0 5 0
0 3 0 1
0 6 0 1
1
end_operator
begin_operator
unstack d b
0
3
0 1 5 0
0 3 0 1
0 6 0 1
1
end_operator
begin_operator
unstack d c
0
3
0 2 5 0
0 3 0 1
0 6 0 1
1
end_operator
begin_operator
unstack d e
0
3
0 3 0 1
0 4 5 0
0 6 0 1
1
end_operator
begin_operator
unstack d f
0
3
0 3 0 1
0 5 5 0
0 6 0 1
1
end_operator
begin_operator
unstack e a
0
3
0 0 6 0
0 4 0 1
0 6 0 1
1
end_operator
begin_operator
unstack e b
0
3
0 1 6 0
0 4 0 1
0 6 0 1
1
end_operator
begin_operator
unstack e c
0
3
0 2 6 0
0 4 0 1
0 6 0 1
1
end_operator
begin_operator
unstack e d
0
3
0 3 6 0
0 4 0 1
0 6 0 1
1
end_operator
begin_operator
unstack e f
0
3
0 4 0 1
0 5 6 0
0 6 0 1
1
end_operator
begin_operator
unstack f a
0
3
0 0 7 0
0 5 0 1
0 6 0 1
1
end_operator
begin_operator
unstack f b
0
3
0 1 7 0
0 5 0 1
0 6 0 1
1
end_operator
begin_operator
unstack f c
0
3
0 2 7 0
0 5 0 1
0 6 0 1
1
end_operator
begin_operator
unstack f d
0
3
0 3 7 0
0 5 0 1
0 6 0 1
1
end_operator
begin_operator
unstack f e
0
3
0 4 7 0
0 5 0 1
0 6 0 1
1
end_operator
0
