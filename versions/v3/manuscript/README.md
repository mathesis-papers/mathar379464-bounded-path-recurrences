# Mathar's path recurrences and alternative proofs of Prodinger's enumeration

Haobo Ma (ChronoAI Pte Ltd, Singapore), followed by Wenlin Zhang (National
University of Singapore; ChronoAI Pte Ltd, Singapore). CC BY 4.0.

## Contents and build

`note.tex` contains the complete note and its bibliography. Compile with
Tectonic 0.17.0:

```sh
tectonic --untrusted note.tex
```

Alternatively run `pdflatex note.tex` twice. No BibTeX step or external figure
is needed. The standard packages are geometry, fontenc, lmodern, amsmath,
amssymb, amsthm, booktabs, array, longtable and hyperref.

`metadata.json` contains the title, abstract, byline, license and source
coverage. `verify.py` and `checks.json` give exact arithmetic checks:

```sh
python3 -m pip install -r requirements.txt
python3 verify.py > checks-new.json
```

The output should equal `checks.json` as a JSON object.

## Mathematical coverage

The written proofs establish Mathar's eight-term recurrence for floor -2
paths for all n >= 7, his six-term recurrence for floor -1 paths for all
n >= 5, the zero-floor five-term recurrence, and the recorded asymptotics.
The Raney enumeration is credited to the earlier literature, including
Prodinger's 2019 preprint and its 2020 revision, and Drube's 2022 treatment.
A polynomial correction gives the generating function for every fixed floor.
The cyclic-cut proof is also given by Asinowski, Hackl and Selkirk (2020/2022,
Proposition 3); the zero-floor sum appears in Barry (2015/2016, Example 16),
and the ordered-tree class follows from Banderier and Flajolet (2002,
Example 3). The note gives full bibliographic locators for these prior results.

Companion corrections are explicit: A127902's printed offset requires N-1
in its binomial sum; its tree interpretation uses ordered rooted trees
without vertex labels; A069271's inverse-power coefficient is 4*(-1)^n,
and one PARI expression needs its extra factor 2 removed. The divided power
formula excludes m = -2, whose value the polynomial expression supplies.
For the zero-free floor -2 sibling A006632, the note proves the recorded
generating function V=t*C^3=t/(1-V)^3, recurrence for N >= 2, and
zero-offset inverse-power characterization.

The finite comparisons cover the three original floors through length 240,
all floor depths 0 through 8 through length 120, powers with multipliers
-5 through 5 and degrees 1 through 16, and the production matrix through
power 30. The differential-operator identities are symbolic exact identities.
These ranges supplement, rather than replace, the all-index written proofs.

## Formal statements

The pinned declaration `D5.S1.Recurrence.MatharA379464.Result.result`
states the eight-term recurrence for the literal floor -2 path count and
all indices n >= 7:

https://github.com/the-omega-institute/trureturing/blob/03150b25ffad5ced1bebdca2efc8c05b484c1eca/D5/S1/Recurrence/MatharA379464/Result.lean

The supporting declarations `three_positive_rotations` and `binomial_counting`
state, respectively, the three-positive-cuts theorem for lists of integers at
most one with sum three, and the floor -2 binomial insertion formula:

https://github.com/the-omega-institute/trureturing/blob/03150b25ffad5ced1bebdca2efc8c05b484c1eca/D5/S1/Words/LatticePaths/MatharA379464Counting.lean

The other floors, analytic asymptotics, all-floor polynomial correction and
companion corrections have the written proofs in the note.
