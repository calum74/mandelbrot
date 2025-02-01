Next idea: When is it reasonable to drop quadratic terms completely?

Recall the basic recurrence relation to calculate the distance $\epsilon_i$ from a reference orbit at iteration $i$:

1. $\epsilon_{i+1} = 2z_i\epsilon_i + \epsilon_i^2 + \delta$

If $\epsilon_i$ is "small" (less than $10^{-16}$ or so), then we can drop the $\epsilon_i^2$ term entirely, since it cannot affect the calculation. This gives

2. $\epsilon_{i+1} = 2z_i\epsilon_i + \delta$

We can then formulate a "jump forward" $n$ steps using linear equations

3. $\epsilon_0 = 0$

4. $\epsilon_{i+n} = A_{i,n}\epsilon_{i} + B_{i,n}\delta$

5. $A_{i, 0} = 1$

6. $B_{i,0} = 0$

**Theorem 1:** (Jump forward any $n$ steps from iteration $i$)

7. &nbsp;&nbsp;&nbsp;&nbsp;$A_{i,n+1} = 2z_{i+n}A_{i,n}$

8. &nbsp;&nbsp;&nbsp;&nbsp;$B_{i,n+1} = 2z_{i+n}B_{i,n}+1$

Proof:

9. &nbsp;&nbsp;&nbsp;&nbsp;$\epsilon_{i+n+1} = A_{i,n+1}\epsilon_{i} + B_{i,n+1}$ from 4

10. &nbsp;&nbsp;&nbsp;&nbsp;$ = 2z_{i+n}\epsilon_{i+n} + \delta$ from 2

11. &nbsp;&nbsp;&nbsp;&nbsp;$= 2z_{i+n}(A_{i,n}\epsilon_i + B_{i,n}\delta)+\delta$ from 4.

12. &nbsp;&nbsp;&nbsp;&nbsp;$= (2z_{i+n}A_{i,n})\epsilon_i + (2z_{i+n}B_{i,n}+1)\delta$

Equate terms in 9. and 12. to give 7. and 8.
$\square$

**Theorem 2:** (Jump forward $n$ steps from any iteration $i$)

13. &nbsp;&nbsp;&nbsp;&nbsp;$A_{i+1,n} = 2z_{i+n}A_{i,n}$

14. &nbsp;&nbsp;&nbsp;&nbsp;$B_{i+1,n} = 2z_{i+n}B_{i,n}+1$

Proof:

15. &nbsp;&nbsp;&nbsp;&nbsp;$\epsilon_{i+n+1} = A_{i+1,n}\epsilon_{i} + B_{i+1,n}$ from 4

Equate terms in 12 and 15. to give 13 and 14.
$\square$

To estimate the size, we need to ensure $\epsilon_i \lt 2^{-52}$

&nbsp;&nbsp;&nbsp;&nbsp;$B_{i,n}\delta \lt 2^{-53}$

&nbsp;&nbsp;&nbsp;&nbsp;$B_{i,n} \lt 2^{-53}/\delta$

In practical terms, we can store 3 extra numbers per iteration in a reference orbit:

* The number of steps $n$ we can jump (may be 0)
* Terms $A_{i,n}$ and $B_{i,n}$ encoding the linear step.

*However* what's the guarantee that we'll ever find an $n>0$ if the series approximation already failed?

## Higher powers

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

