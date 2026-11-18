# Bounded-suboptimal, anytime and incremental search

The searches that give up the optimality of the path for a bound on its
cost, those that improve a first path until they prove it optimal, and those
that repair a search after some costs have changed are treated in [chapter 6](https://dmeoli.github.io/PBeFS/manual/chapters/06-bounded-suboptimal.html) of
the [manual](https://dmeoli.github.io/PBeFS/manual/index.html): the weight $w$ and the rules that keep the bound
$w C^*$ in any order of expansion, the duplicate dropping and the heuristics
for which it keeps the bound, in [section 6.2](https://dmeoli.github.io/PBeFS/manual/chapters/06-bounded-suboptimal.html#sec-6-2); weighted A* and weighted IDA* in
[section 6.3](https://dmeoli.github.io/PBeFS/manual/chapters/06-bounded-suboptimal.html#sec-6-3); the weighted and anytime parallel searches in [section 6.4](https://dmeoli.github.io/PBeFS/manual/chapters/06-bounded-suboptimal.html#sec-6-4); AWA*, ARA*
and the Anytime Restarting PBNF in [section 6.5](https://dmeoli.github.io/PBeFS/manual/chapters/06-bounded-suboptimal.html#sec-6-5); the optimistic search, the
Improved Optimistic Search, the focal searches and the Dynamic Potential Search
in [section 6.6](https://dmeoli.github.io/PBeFS/manual/chapters/06-bounded-suboptimal.html#sec-6-6); the batch weighted A* in [section 6.7](https://dmeoli.github.io/PBeFS/manual/chapters/06-bounded-suboptimal.html#sec-6-7); LPA* and D* Lite with the
changes of the costs told as events in [section 6.8](https://dmeoli.github.io/PBeFS/manual/chapters/06-bounded-suboptimal.html#sec-6-8); and LPA* on the workers of
HDA*, with its hypotheses, the theorems that every episode ends with an
optimal path whatever the interleaving, the window on the least key and the
choice between repair and restart, in [section 6.9](https://dmeoli.github.io/PBeFS/manual/chapters/06-bounded-suboptimal.html#sec-6-9).
