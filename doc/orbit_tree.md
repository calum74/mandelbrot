# Orbit tree algorithm

Invented by Calum Grant
January 15-18 2025
Citation: github.com/calum74/mandelbrot/doc/orbit_tree.md

## Overview

*Orbit-tree* is a novel and efficient algorithm for rendering the Mandelbrot set. It builds on the work of perturbation orbits and Taylor series.

As its name suggests, it adds the capability to efficiently *branch* an orbit, so when we notice that an orbit becomes unstable, we can reset the orbit and branch the orbit into 4 nearby orbits. This is efficient because the new orbits only need to be calculated from the current iteration, and not iteration 0. This gives the capability to create reference orbits arbitrarily close to any point in the image.

This data structure is called an *orbit-tree* because the orbits form a tree-like structure. The "length" of a branch is the number of iterations in the orbit before it becomes unstable or escapes. The "height" of the tree is the number of iterations. The "root" of the tree is the reference orbit at the center of the image. A "leaf" of the tree represents a small region of the image.

Orbit-trees contain Taylor series terms for the evaluation of an epsilon (dz) relative to the current branch. This Taylor series is used to quickly compute the distance of a point from the reference orbit for any iteration.

To evaluate a point using the orbit-tree, we can look down the tree into in each branch, noting that the number of branches is fixed at O(ln W). (W is the width of the image in pixels.) Using a binary search, we can locate the escape iteration in a branch in O(ln N) time. (N is the number of iterations in the branch.)  Furthermore, the escape iteration tends to be towards the top of the tree, and we can use the previous foud escape iteration as a starting point.

This algorithm solves the fundamental problem with the simple Taylor series technique, which either creates too many reference orbits (expensive), or too few reference orbits resulting in series divergence and too few skip iterations.

This algorithm offers an improvement over traditional Taylor series techniques, because the Taylor series has a tendency to diverge, particularly for complex images. So the number of skip-iterations can still be low relative to the depth of the image.

## Branching criteria

Ideally we want to allow branches to become as long as possible

Note that the Taylor series terms can become very large, so a bespoke number implementation is needed to handle numbers of very large exponents. Without this, the length of each branch is limited to 100 or so.

## How to branch

## Outline

## Implementation notes

To implement orbit-tree, you will need

* High-precision numbers
* High-exponent numbers
* Complex numbers
