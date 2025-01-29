# The orbit-tree algorithm for rendering the Mandelbrot set

Original work by Calum Grant, January 2025

Cite as https://github.com/calum74/mandelbrot

## Overview

*Orbit-tree* is a novel algorithm for rendering the Mandelbrot set, that builds upon perturbation theory and Taylor series expansion. Recall that we can calculate terms for $\epsilon_i = A_i\delta + B_i\delta^2+C_i\delta^3 + ...$ that allow us to approximate the distance of any orbit in the Mandelbrot set at a distance $\delta$ from a *reference orbit*. This provides a significant speedup since it avoids iterating each point.

The main problem with the Taylor series expansion is that if the reference orbit is too far away from the point being calculated, then the terms of the Taylor series diverge. This is particularly a problem for very complex images that have a very large range of iterations.

The orbit-tree algorithm makes it efficient to create reference orbits arbitarily close to any point in the image, which means that the Taylor series is always optimal. 

The innovation is that we can branch a reference orbit as soon as it loses precision for its region. This creates 4 sub-orbits of higher precision, and we carry on branching reference orbits until we reach a minimum size.
The contribution is that we can branch not just the relative orbit, but also the terms in the Taylor series, in O(1) time without needing to recompute the entire series.

Once constructed, the orbit-tree ensures that all Taylor series are valid for all points in the region.

## Mathematical foundations

**Definition 1:** Let $z,z',\Delta \in \mathbb{C}$ be reference orbits, such that $z'=z+\Delta$, and $z_i$ and $z'_i$ are the $i^{th}$ iteration of $z$ and $z'$ of their respective orbits. 

Here we mean an *orbit* to be the sequence of numbers defined by $z_{i+1}=z_i^n+z$, or similar, used to compute the Mandelbrot set.

**Definition 2:** Let $\delta \in \mathbb{C}$ be the distance of a point $z''$ from $z$, and $\delta' \in \mathbb{C}$ be the distance of $z''$ from $z'$. 

It follows that $z+\delta=z'+\delta'$, therefore $z+\delta=z+\Delta+\delta'$, therefore $\delta=\Delta+\delta'$.

**Definition 3:** $\epsilon_i \in \mathbb{C} = z''_i - z_i$, and $\epsilon'_i \in \mathbb{C} = z''_i - z'_i$ 

$\epsilon_i$ represents the distance of orbit $z''$ from the reference orbit $z$, and $\epsilon'_i$ is the distance of orbit $z''$ from $z'$.

**Definition 4:** Let $\Epsilon_i \in \mathbb{C}$ be the distance between the reference orbits, so $\Epsilon_i = \epsilon_i - \epsilon'_i$

**Definition 5:** For terms $A_i,B_i,C_i,D_i \in \mathbb{C}$, $\epsilon_i = A_i\delta + B_i\delta^2 + C_i\delta^3 + D_i\delta^4 + O(\delta^5)$

**Definition 6:** For terms $A'_i,B'_i,C'_i,D'_i \in \mathbb{C}$, $\epsilon'_i = A'_i\delta' + B'_i\delta'^2 + C'_i\delta'^3 + D'_i\delta'^4 + O(\delta'^5)$

These are the familiar Taylor series expansions for $\epsilon_i, \epsilon'_i$ which we can use to skip iterations when computing the escape iteration in a Mandelbrot set.

**Theorem 1:** (Translating a Taylor series by $\Delta$)

$A'_i = A_i + 2B_i\Delta + 3C_i\Delta^2 + 4D_i\Delta^3$

$B'_i = B_i + 3C_i\Delta + 6D_i\Delta^2$

$C'_i = C_i + 4D_i\Delta$

$D'_i = D_i$

*Proof:*

$\Epsilon_i = A_i\Delta + B_i\Delta^2 + C_i\Delta^3 + D_i\Delta^4 + O(\Delta^5)$

$\epsilon'_i = A'_i\delta' + B'_i\delta'^2 + C'_i\delta'^3 + D'_i\delta'^4 + O(\delta'^5)$ (from Definition 6)

$= \epsilon_i - \Epsilon_i$ (from Definition 4)

$=  A_i\delta + B_i\delta^2 + C_i\delta^3 + D_i\delta^4 + O(\delta^5) - (A_i\Delta + B_i\Delta^2 + C_i\Delta^3 + D'_i\Delta^4 + O(\Delta^5))$

$= A_i(\delta'+\Delta) + B_i(\delta'+\Delta)^2 + C_i(\delta'+\Delta)^3 + D_i(\delta'+\Delta)^4 + O((\delta'+\Delta)^5) - (A_i\Delta + B_i\Delta^2 + C_i\Delta^3 + D_i\Delta^4 + O(\Delta^5))$

$= A_i\delta'+A_i\Delta + B_i(\delta'^2 + 2\delta'\Delta + \Delta^2) + C_i(\delta'^3+3\delta^2\Delta + 3\delta'\Delta^2 + \Delta^3) + D_i(\delta'^4+4\delta'\Delta^3 + 6\delta'^2\Delta^2 + 4\delta'\Delta^3 + \Delta^4) + O((\delta'+\Delta)^5) - (A_i\Delta + B_i\Delta^2 + C_i\Delta^3 + D_i\Delta^4 + O(\Delta^5))$

$= A_i\delta'+A_i\Delta + B_i\delta'^2 + 2B_i\delta\Delta + B_i\Delta^2 + C_i\delta^3+3\delta^2\Delta + 3C_i\delta'\Delta^2 + C_i\Delta^3 + D_i\delta'^4+4D_i\delta'\Delta^3 + 6D_i\delta'^2\Delta^2 + 4D_i\delta'\Delta^3 + D_i\Delta^4 + O((\delta'+\Delta)^5) - A_i\Delta - B_i\Delta^2 - C_i\Delta^3 - D_i\Delta^4 - O(\Delta^5)$

$= A_i\delta'+A_i\Delta + B_i\delta'^2 + 2B_i\delta'\Delta + B_i\Delta^2 + C_i\delta^3+3\delta'^2\Delta + 3C_i\delta\Delta^2 + C_i\Delta^3 + D_i\delta'^4+4D_i\delta'\Delta^3 + 6D_i\delta'^2\Delta^2 + 4D_i\delta'\Delta^3 + D_i\Delta^4 + O((\delta'+\Delta)^5) - A_i\Delta - B_i\Delta^2 - C_i\Delta^3 - D_i\Delta^4 - O(\Delta^5)$

$= (A_i+ 2B_i\Delta + 3C_i\Delta^2+ 4D_i\Delta^3)\delta' + (B_i + 3C_i\Delta+ 6D_i\Delta^2)\delta'^2  + (C_i +4D_i\Delta)\delta'^3 + D_i\delta'^4 + O((\delta'+\Delta)^5)$

Equating terms in $\delta'$, $\delta'^2$, $\delta'^3$ and $\delta'^4$, we get Theorem 1. $\square$

What this means is that we can translate any Taylor series expansion of $\epsilon_i$ to a new orbit without recomputing the entire series. If we have $A_i,B_i,C_i,D_i$ from one reference orbit, we can compute $A'_i,B'_i,C'_i,D'_i$ for a second reference orbit at a distance $\Delta$.

We can generalise Theorem 1 to any power of $\delta$, and the terms are from Pascal's triangle. (Exercise for the reader.)

**Theorem 2:** (Translating a relative orbit by $\Delta$ and $\Epsilon_i$)

For a relative orbit $z' = z + \delta, z'_i = z_i + \epsilon_i$, an orbit $z'' = z'+\Delta, z''_i = z'_i+\Epsilon_i$, then $z'' = z+(\delta+\Delta), z''_i = z_i + (\epsilon_i + \Epsilon_i)$.

Proof: $z'' = z'+\Delta = (z + \delta)+\Delta = z+(\delta+\Delta)$. $z''_i = z'_i+\Epsilon_i = (z_i + \epsilon_i) + \Epsilon_i = z_i + (\epsilon_i + \Epsilon_i)$. $\square$

This just says that relative orbits can be easily translated simply by adding $\Delta$ and $\Epsilon_i$.

Theorem 1 and Theorem 2 together allow us to translate a Taylor series reference orbit by $\Delta$ and $\Epsilon_i$, assuming we can accurately calculate $\Epsilon_i$, which we do by checking that the terms in the Taylor series don't diverge too much.

## Algorithm outline

### Orbit-trees

An orbit-tree consists of a *root*, *branches*, and *leaves*. The root is an orbit at the center of the image, starting at iteration 0.

All branches have one parent and 4 children, except for the root which has no parent, and leaves which have no children.

Branches represent a range of iterations, where the last iteration of its parent is the first iteration of the child.

Branches store $z'_i$ and terms $A_i,B_i,C_i,D_i$ for each iteration $i$ in an array. The range of a branch represents the iterations for which the Taylor series is valid within its radius $r \in \mathbb{R}$. A branch has a $\Delta$ to its parent branch.

### Creating the orbit-tree

Starting at the root, compute each iteration in the branch, using the well-known equations, see [Wikipedia](https://en.wikipedia.org/wiki/Plotting_algorithms_for_the_Mandelbrot_set).

We iterate until:

* The terms $A_i, B_i, C_i, D_i$ diverge for radius $r$
* $i$ exceeds the iteration limit
* $z'_i$ escapes

We will then either:

1. Create child branches, and repeat
2. Compute all pixels in a leaf branch

### Creating child branches

Compute deltas $\Delta_a, \Delta_b,\Delta_c,\Delta_d$ representing the relative delta to the center of each branch, of length $r/2$.

For each branch, compute the new terms $A_i, B_i, C_i, D_i, \epsilon_i$ using the formulae

$A'_i = A_i + 2B_i\Delta + 3C_i\Delta^2 + 4D_i\Delta^3$

$B'_i = B_i + 3C_i\Delta + 6D_i\Delta^2$

$C'_i = C_i + 4D_i\Delta$

$D'_i = D_i$

$\epsilon'_i = \epsilon_i + A_i\Delta + B_i\Delta^2 + C_i\Delta^3 + D_i\Delta^4$

(If using 3 terms, $D_i$ will be 0).

### Computing the escape iteration

Compute $\delta'$ of each pixel $z''$ to its branch.

We can use a binary search in the branch to see if the orbit has escaped at that iteration. Evaluate the orbit value $z''_i = z'_i + A_i\delta' + B_i\delta'^2 + C_i\delta'^3 + D_i\delta'^4$ using the terms from the branch and test if it has escaped.

If a branch is fully escaped at $\delta'$, we need to check the *parent* branch, which involves subtracting the relevant $\Delta$ from $\delta'$ each time we transition branches. Since the height of the orbit-tree is $O(ln(n))$, it should take at most $O(ln(n))$ steps to calculate each point, and we can use the previous iteration count as a starting point.

If may also be that a pixel is not escaped even at the leaf branch. In that case, we need to carry on iterating the pixel $z''_i$ until it escapes or we reach the iteration limit.

## Performance

To do!

