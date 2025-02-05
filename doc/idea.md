*This is a work in progress*

Next idea: When is it reasonable to drop quadratic terms completely?

# Bilinear approximation

Recall the basic recurrence relation to calculate the distance $\epsilon_i$ from a reference orbit at iteration $i$:

1. $\epsilon_{i+1} = 2z_i\epsilon_i + \epsilon_i^2 + \delta$

If $\epsilon_i$ is "small" (less than $10^{-16}$ or so), then we can drop the $\epsilon_i^2$ term entirely, since it cannot affect the calculation. This gives

2. $\epsilon_{i+1} = 2z_i\epsilon_i + \delta$

We can then formulate a "jump forward" $n$ steps from iteration $i$ using linear equations

3. $\epsilon_0 = 0$
4. $\epsilon_{i+n} = A_{i,n}\epsilon_{i} + B_{i,n}\delta$
5. $A_{i, 0} = 1$
6. $B_{i,0} = 0$

**Theorem 1:**

7. &nbsp;&nbsp;&nbsp;&nbsp;$A_{i,n+1} = 2z_{i+n}A_{i,n}$
8. &nbsp;&nbsp;&nbsp;&nbsp;$B_{i,n+1} = 2z_{i+n}B_{i,n}+1$

Proof:

9. &nbsp;&nbsp;&nbsp;&nbsp;$\epsilon_{i+n+1} = A_{i,n+1}\epsilon_{i} + B_{i,n+1}$ from 4
10. &nbsp;&nbsp;&nbsp;&nbsp;$ = 2z_{i+n}\epsilon_{i+n} + \delta$ from 2
11. &nbsp;&nbsp;&nbsp;&nbsp;$= 2z_{i+n}(A_{i,n}\epsilon_i + B_{i,n}\delta)+\delta$ from 4.
12. &nbsp;&nbsp;&nbsp;&nbsp;$= (2z_{i+n}A_{i,n})\epsilon_i + (2z_{i+n}B_{i,n}+1)\delta$

Equate terms in 9. and 12. to give 7. and 8.
$\square$

**Theorem 2:** (Composition)

13. $A_{i,n+m} = A_{i,n}A_{i+n,m}$
14. $B_{i,n+m} = A_{i+n,m}B_{i,n}+B_{i+n,m}$

Proof:

15. $\epsilon_{i+n+m} = A_{i+n,m}\epsilon_{i+n} + B_{i+n,m}\delta$
16. $= A_{i,n+m}\epsilon_i + B_{i,n+m}\delta$
17. $= A_{i+n,m}(A_{i,n}\epsilon_i + B_{i,n}\delta) + B_{i+n,m}\delta$ expanding 4 into 15
18. $= A_{i,n}A_{i+n,m}\epsilon_i + (A_{i+n,m}B_{i,n}+B_{i+n,m})\delta$ expanding 17

Equating terms in 16 and 18, we get 13 and 14. $\square$

**Theorem 3:**

21. $A_{i+1,m} = z_{i+m}A_{i,m}/z_i$
22. $B_{i+1,m} = 2z_{i+m}B_{i,m}+1 - A_{i+1,m}$

Proof:

From 7,

23. $A_{i,1} = 2z_i$

Using Theorem 2 with n=1,

24. $A_{i,m+1} = 2z_{i+m}A_{i,m}$

25. $A_{i+1,m} = A_{i,m+1}/A_{i,n} = 2z_{i+m}A_{i,m}/A_{i,1} = 2z_{i+m}A_{i,m}/2z_i = z_{i+m}A_{i,m}/z_i$

Using Theorem 2 with n=1

26. $B_{i,m+1} = A_{i+1,m}B_{i,1}+B_{i+1,m} = A_{i+1,m} + B_{i+1,m}$

27. $B_{i+1,m} = B_{i,m+1} - A_{i+1,m}$

$\square$

Corollary 1 gives us a way to calculate a way to jump forward $m$ steps for each iteration.


## Precision of bilinear approximation

Define $\epsilon \in R$ to be the smallest number representable using a floating point number such that if $b<\epsilon a$ then $a \equiv a+b$.

Note we are using $\epsilon$ to indicate floating point precision, whereas $\epsilon_i$ is the distance from the reference orbit. Confusing!

**Theorem 4:**

If

$|A| < \epsilon |z_i||\epsilon_j|^{-1}$

$|A| < \epsilon^{1/2}|\delta|^{1/2}|\epsilon_j|^{-1}/2$

$|B| < \epsilon |z_i||\delta|^{-1}$

$|B|<\epsilon^{1/2}|\delta|^{-1/2}/2$ 

then the linearization is valid.

Proof:

The linearization remains valid so long as

$2z_i\epsilon_i + \epsilon_i^2 + \delta \equiv 2z_i\epsilon_i + \delta$

$\impliedby |\epsilon_i^2| < \epsilon|2z_i\epsilon_i + \delta| < 2\epsilon|z_i\epsilon_i| + \epsilon|\delta| = 2\epsilon|z_i||\epsilon_i| + \epsilon|\delta|$ by the $\triangle$ rule

$\impliedby |\epsilon_i|^2 < 2\epsilon |z_i||\epsilon_i|$ and $|\epsilon_i|^2 < \epsilon|\delta|$ Fishy!!!

$\iff |\epsilon_i| < 2\epsilon |z_i|$ and $|\epsilon_i| < \epsilon^{1/2}|\delta|^{1/2}$

But $\epsilon_i = A\epsilon_j + B\delta$ so we need to derive limits for $A$ and $B$ in terms of $\epsilon_j$, $z_i$ and $\delta$.

$\iff |A\epsilon_j + B\delta| < 2\epsilon |z_i|$ and $|A\epsilon_j + B\delta| < \epsilon^{1/2}|\delta|^{1/2}$

$\impliedby |A\epsilon_j| + |B\delta| < 2\epsilon |z_i|$ and $|A\epsilon_j| + |B\delta| < \epsilon^{1/2}|\delta|^{1/2}$ 

$\impliedby |A\epsilon_j| < \epsilon |z_i|$ and $|B\delta| < \epsilon |z|$ and $|A\epsilon_j| < \epsilon^{1/2}|\delta|^{1/2}/2$ and $|B\delta|<\epsilon^{1/2}|\delta|^{1/2}/2$ 

$\iff |A||\epsilon_j| < \epsilon |z_i|$ and $|B||\delta| < \epsilon |z_i|$ and $|A||\epsilon_j| < \epsilon^{1/2}|\delta|^{1/2}$ and $|B||\delta|<\epsilon^{1/2}|\delta|^{1/2}$ 

$\iff |A| < \epsilon |z_i||\epsilon_j|^{-1}$ and $|B| < \epsilon |z_i||\delta|^{-1}$ and $|A| < \epsilon^{1/2}|\delta|^{1/2}|\epsilon_j|^{-1}/2$ and $|B|<\epsilon^{1/2}|\delta|^{-1/2}/2$  

$\square$

The question is, how long can the series go on for? For this, we'll need an idea of the maximum value of $|\epsilon_j|$ given by

$|\epsilon_j| < |A|max(|\epsilon_k|) + |B|max(|\delta|)$ 

where $max|\epsilon_k|$ is the maximum value of $|\epsilon_k|$ from the previous branch.

This gives the maximum permitted value for $|\epsilon_j|$ which in turn gives a maximum permitted value for the coefficients $|A|$ and $|B|$ at each iteration $i$. We iterate the branch until the terms $|A|$ and $|B|$ reach their maximum values permitted by Theorem 4, and this is the length of the branch.

Notes:


Theorem:

$|A|^2 < \epsilon^2|z|^2|\epsilon_j|^{-2}$ and $|B|^2 < \epsilon^2|z|^2|\delta|^{-2}$

Proof:
$|\epsilon_i| < \epsilon|z|$ approximately, so just watch out for $z$ close to 0.

$|A\epsilon_j + B\delta| < \epsilon|z|$ so approximately

$|A||\epsilon_j| < \epsilon|z|$ and $|B||\delta| < \epsilon|z|$

$|A| < \epsilon|z||\epsilon_j|^{-1}$ and $|B| < \epsilon|z||\delta|^{-1}$

$\square$

## Quadratic approximation

What if we instead created a quadratic approximation for skip-forward:

$\epsilon_{i+n} = A_{i,n}\epsilon_i + A'_{i,n}\epsilon_i^2 + B_{i,n}\delta + B'_{i,n}\delta^2 + C_{i,n}\epsilon\delta + O(\epsilon_i^3)$

Where $A_{i,n}$ is a coefficient for skipping $n$ steps forward from iteration $i$.

$A_{i,0} = 1$, $A'_{i,0} = B_{i,0} = B'_{i,0} = C_{i,0} = 0$

**Theorem 5** (Computing quadratic terms)

$A_{i,n+1} = 2z_{i+n}A_{i+n}$

$A'_{i,n+1} = 2z_{i+n}A'_{i+n} + A_{i+n}^2$

$B_{i,n+1} = 2z_{i+n}B_{i,n}+1$

$B'_{i,n+1} = 2z_{i+n}B'_{i,n} + B_{i,n}^2$

$C_{i,n+1} = 2z_{i,n}C_{i,n}+2A_{i,n}B_{i,n}$

Proof:

1. $\epsilon_{i+n+1} = A_{i,n+1}\epsilon_i + A'_{i,n+1}\epsilon_i^2 + B_{i,n+1}\delta + B'_{i,n+1}\delta^2 + C_{i,n+1}\epsilon\delta + O(\epsilon_i^3)$

from the definition. We also have

2. $\epsilon_{i+n+1} = 2z_{i+n}\epsilon_{i+n} + \epsilon_{i+n}^2 + \delta$

3. $= 2z_{i+n}(A_{i,n}\epsilon_i + A'_{i,n}\epsilon_i^2 + B_{i,n}\delta + B'_{i,n}\delta^2 + C_{i,n}\epsilon\delta) + (A_{i,n}\epsilon_i + A'_{i,n}\epsilon_i^2 + B_{i,n}\delta + B'_{i,n}\delta^2 + C_{i,n}\epsilon\delta)^2 + \delta$

4. $= 2z_{i+n}(A_{i,n}\epsilon_i + A'_{i,n}\epsilon_i^2 + B_{i,n}\delta + B'_{i,n}\delta^2 + C_{i,n}\epsilon\delta) + (A_{i,n}\epsilon_i + B_{i,n}\delta)^2 + \delta + O(\epsilon^3)$ 

5. $= (2z_{i+n}A_{i+n})\epsilon_i + (2z_{i+n}A'_{i+n} + A_{i+n}^2)\epsilon_i^2 + (2z_{i+n}B_{i,n}+1)\delta + (2z_{i+n}B'_{i,n} + B_{i,n}^2)\delta^2 + (2z_{i,n}C_{i,n}+2A_{i,n}B_{i,n})\epsilon\delta + O(\epsilon^3)$

We can then equate terms in 1 and 5.

$\square$

We can also look at the size of the "residual" from equation 4 above.

This approximation is valid provided that $|\epsilon_i|^2 < \epsilon$, which should give a longer branch. It is unclear whether the additional effort in computing extra terms is outweighed by having a longer branch.

# Relocating the tree

Imagine that you have a reference orbit $z$, a relative orbit $z'=z+\delta$, and you want to create a third orbit $z''=c'+\delta'$ close to the second orbit.

At iteration $i$, we have $\epsilon_i$ is the distance from the relative orbit to the reference orbit. We also have $\epsilon'_i$ as the distance between the second orbit $z'$ and the third orbit $z''$. It follows that the distance between $z''$ and $z$ is $\epsilon_i + \epsilon'_i$.

We also want to compute terms relative to the second orbit $z'$.



# Turning this idea into an algorithm

The general idea is that we can split the entire iteration into a series of linear big-steps. We can build this on top of perturbation and series-approximation, so for any orbit where we notice we are iterating forwards, we can "detect" that epsilon is small and use this to accelerate iteration.

Another idea is to implement some kind of "rolling hash" so we compute a "skip forward 5" term for every entry. When can we apply then?

## Constructing the tree

Start by computing the orbit from iteration 0, Calculating $A_{0,n}$ and $B_{0,n}$ as far as the precision of $\epsilon_i$ allows.

Then we'll create 4 further orbits, where each sub-orbit (at a distance $\Delta$ from the original) can be started using $\Epsilon_i$, $z'_i$ from the reference orbit. Again we'll compute $A_{n_1,n_2}, B_{n_1,n_2}$ for the sub-orbit, which can be done precisely. Since the maximum $\delta$ of the sub-orbit is less than the original orbit, we can calculate this precisely.

Each branch contains:

* $I \in N$, its base iteration
* $\delta_b$, the distance to the high precision reference orbit
* $\epsilon_I$ the starting epsilon for the branch, relative to the reference orbit, constructed by jumping forward from the previous branches
* $J \in N$, the maximum iteration we can skip forward to
* $A_{i}, B_{i} \in C$ for each iteration, to skip forward to iteration $i$.
* $z_{i} \in C$ for each iteration, the orbital value at each iteration $i$, computed iteratively from $\delta_b$ and $\epsilon_I$.

## Calculating escape iteration

To calculate $z_i$ for a point $p$, we calculate its $\delta$ relative to the first reference orbit. Then we jump forward $J$ places to the second orbit, which maintains our $\epsilon$.

This is very similar to when we branch an orbit. We need to calculate its $\epsilon$ from iteration 0.

