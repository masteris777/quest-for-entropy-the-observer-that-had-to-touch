# The Observer That Had To Touch — companion code

Everything the article claims, runnable from scratch.

## Run it

```
pip install numpy scipy sympy
python run_all.py
```

`run_all.py` runs every experiment and checks each claim the article makes,
naming the sentence each check is testing. The two fits are the slow part;
`python run_all.py --quick` skips them and runs the proof alone.

## What is in here

The article has two halves — the measuring that failed, and the proof that it
had to. The code follows the same split.

**The measuring**

| file | what it does |
|---|---|
| `engine.py` | the deterministic machine everything here runs on |
| `single_moment_fit.py` | fits the quantum probability rule to one snapshot — of the machine, of a chaotic system, and of pure random numbers |
| `two_moment_quantum.py` | the harder question: two snapshots, counted frequencies against the quantum collapse formula |
| `the_rotor.py` | a plain wheel turning at a steady rate, fitted to the same twenty comparisons |

The point of `single_moment_fit.py` is that it passes on all three inputs.
A test a random number generator passes is not a test, and finding that out
is what sent the rest of this episode looking for a proof instead.

**The proof**

| file | what it does |
|---|---|
| `the_proof.py` | five lemmas in exact rational and symbolic arithmetic — no floating point anywhere in it |

Its lemmas: the joint distribution exists constructively; every one of the
eight possible outcome-triples obeys the bound (enumerated, not argued);
pointwise implies on-average; the quantum instance in closed form; and — the
one the article ends on — an exhaustive search over *every* passive setting of
the same machine, none of which reproduce what the touching version does.

`article.md` is the piece itself. The interactive page it links is separate and
needs no code from here.

## Scope

The article's "What this does NOT claim" section is the scope fence, and it
ships with the code in `article.md`. Short version: this is a no-go result
about a class of approaches, plus an honest account of a route closing. It
derives nothing.

## Licence

MIT for the code. The article text is © Marijus Masteika.
