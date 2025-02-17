# BLA2

## Introduction

The BLA2 algorithm improves on existing bivariate linear approximation methods by computing higher order terms, and by utilising a novel(?) technique called *orbit shifting* to always locate an efficient reference orbit for bilinear approximation.

### Notation quick reference

$\epsilon_{nij}$ is the distance between two orbits $i$ and $j$ at iteration $n$.

$\delta_{ij}$ is the distance between two orbits $i$ and $j$ at iteration 1.

$A_{nmi}$ and $B_{nmi}$ are polynomial coefficients used to jump from iteration $n$ to iteration $m$ in the region of orbit $i$.

## Basic definitions

Firstly let's establish some basic notation. When computing a Mandelbrot set, we need to compute whether each point $c_i \in \mathbb{C}$ *escapes*, which means defining a sequence of complex numbers $z_{n,i}$ calculated as follows:

**Definition 1:** Let $z_{ni} \in \mathbb{C}$ be the sequence of points in orbit $i$, defined by

$$\begin{align}z_{0i} &= 0 \newline
z_{n+1,i} &= z_{ni}^2 + c_i\end{align}$$

for some point $c_i \in \mathbb{C}$.
$\square$

Since we will be referring to multiple orbits, we label each orbit $i,j,k$ etc.

**Definition 2:** Let $\epsilon_{nij} \in \mathbb{C} = z_{ni} - z_{nj}$ be the distance between two orbits $i$ and $j$ at iteration $n$. $\square$

**Definition 3:** Let $\delta_{ij} \in \mathbb{C} = c_i - c_j$ to be the distance of two orbits. $\square$

Not to be confused with tensor notation and the Kronecker delta. This notation is based on [Wikipedia notation](https://en.wikipedia.org/wiki/Plotting_algorithms_for_the_Mandelbrot_set).

**Theorem 1:** (Basic properties of $\delta_{ij}$ and $\epsilon_{nij}$)

$$\begin{align}\epsilon_{nij} &= - \epsilon_{nji}\newline
\epsilon_{nii} &= 0\newline
\epsilon_{nik} &= \epsilon_{nij} + \epsilon_{njk} \newline
\delta_{ij} &= \epsilon_{1ij}\newline
\delta_{ii} &= 0\newline
\delta_{ij} &= -\delta_{ji}\newline
\delta_{ik} &= \delta_{ij} + \delta_{jk}
\end{align}$$

Proof: These follow straightforwardly from Definition 2 and 3.

$$\begin{align}
\epsilon_{nij} &= z_{ni} - z_{nj} \newline
&= - (z_{nj} - z_{ni}) \newline 
&= - \epsilon_{nji}\newline
\epsilon_{nii} &= z_{ni} - z_{ni} \newline
&= 0 \newline
\epsilon_{nik} &= z_{ni} - z_{nk} \newline
&=z_{ni} - z_{nj} + z_{nj} - z_{nk} \newline
&= (z_{ni} - z_{nj}) + (z_{nj} - z_{nk}) \newline
&= \epsilon_{nij} + \epsilon_{njk}
\end{align}$$

The $\delta_{ij}$ case follows similiarly, or from $\delta_{ij}=\epsilon_{1ij}$. $\square$

Next we'll show a way to iteratively calculate $\epsilon_{n,i,j}$. Traditional perturbation theory only considers a single reference orbit, but we'll extend this to an arbitrary pair of orbits $i$ and $j$.

**Theorem 2:** (Iterative properties of $\epsilon_{nij}$)

$$\begin{align}\epsilon_{n+1,i,j} = 2z_{nj}\epsilon_{nij}+ \epsilon_{nij}^2 + \delta_{ij}\newline\end{align}$$

Proof:
From Definition 2, 
$$\begin{align}
\epsilon_{n+1,i,j} &= z_{n+1,i} - z_{n+1,j} \newline
\end{align}$$

From Definition 1,

$$\begin{align}
&= (z_{ni}^2 + c_i) - (z_{nj}^2 + c_j) \newline
 &= z_{ni}^2 - z_{nj}^2 + (c_i - c_j) \newline
\end{align}$$

From Definition 3,

$$\begin{align}
 &= z_{ni}^2 - z_{nj}^2 + \delta_{ij} \newline
\end{align}$$

From Definition 2,

$$\begin{align}
 &= (z_{nj} + \epsilon_{nij})^2 - z_{nj}^2 + \delta_{ij} \newline
 &= (z_{nj}^2 + 2z_{nj}\epsilon_{nij}+ \epsilon_{nij}^2) - z_{nj}^2 + \delta_{ij}\newline
 &= 2z_{nj}\epsilon_{nij} + \epsilon_{nij}^2 + \delta_{ij}\end{align}$$

$\square$

Making $i$ and $j$ implicit we get the familiar

$$\begin{align}\epsilon_{n+1} = 2z_{n}\epsilon_{n}+ \epsilon_{n}^2 + \delta\end{align}$$ 

## Polynomial approximation

In this next definition, we'll create a formula for skipping $m-n$ iterations using a quadratic step. Linear approximation (BLA) falls out of the quadratic approximation.

**Definition 5:** (Quadratic approximation)

Let $A_{nmj}$, $A'_{nmj}$, $B_{nmj}$, $B'_{nmj}$, $C_{nmj} \in \mathbb{C}$ be constants such that for orbits $i$ and $j$, the distance between the orbits at iteration $m$ is approximated by

$$\epsilon_{mij} = A_{nmj}\epsilon_{nij} + A'_{nmj}\epsilon_{nij}^2 + B_{nmj}\delta_{ij} + B'_{nmj}\delta_{ij}^2 + C_{nmj}\epsilon_{nij}\delta_{ij} + O(\epsilon_{nij}^3)$$
$\square$

Note that

$$\begin{align}A_{nnj}=1, A'_{nnj} = B_{nnj} = B'_{nnj} = C_{nnj} = 0\end{align}$$

Since

$$\begin{align}\delta_{ij} = O(\epsilon_{nij})\end{align}$$

follows from Theorem 2, we'll gather up all higher order terms and represent them as $O(\epsilon_{nij}^3)$.

**Theorem 3:** (Computation of quadratic terms.)

$$\begin{align}
A_{n,m+1,j} &= 2z_{mj}A_{nmj}\newline

A'_{n,m+1,j} &= 2z_{mj}A'_{nmj} + A_{nmj}^2\newline

B_{n,m+1,j} &= 2z_{mj}B_{nmj}+1\newline

B'_{n,m+1,j} &= 2z_{mj}B'_{nmj} + B_{nmj}^2\newline

C_{n,m+1,j} &= 2z_{mj}C_{nmj}+2A_{nmj}B_{nmj}\end{align}$$

Proof:

$$\begin{align}

\epsilon_{m+1,i,j} &= A_{n,m+1,j}\epsilon_{nij} + A'_{n,m+1,j}\epsilon_{nij}^2 \newline
&+ B_{n,m+1,j}\delta_{ij} + B'_{n,m+1,j}\delta_{ij}^2 \newline 
&+ C_{n,m+1,j}\epsilon_{nij}\delta_{ij} + O(\epsilon_{nij}^3)\newline

\epsilon_{m+1,i,j} &= 2z_{mj}\epsilon_{mij} + \epsilon_{mij}^2 + \delta_{ij}\newline

&= 2z_{mj}(A_{nmj}\epsilon_{nij} + A'_{nmij}\epsilon_{nij}^2 + B_{nmij}\delta_{ij} + B'_{nmj}\delta_{ij}^2 + C_{nmj}\epsilon_{nij}\delta_{ij} + O(\epsilon_{nij}^3)) \newline
&+ (A_{nmij}\epsilon_{nij} + A'_{nmj}\epsilon_{nij}^2 + B_{nmj}\delta_{ij} + B'_{nmj}\delta_{ij}^2 + C_{nmj}\epsilon_{nij}\delta_{ij} + O(\epsilon_{nij}^3))^2 + \delta_{ij}\newline
\end{align}
$$
Gather the higher order terms into $O(\epsilon_{nij}^3)$:
$$\begin{align}
&= 2z_{mj}(A_{nmj}\epsilon_{nij} + A'_{nmj}\epsilon_{nij}^2 + B_{nmj}\delta_{ij} + B'_{nmij}\delta_{ij}^2 + C_{nmj}\epsilon_{nij}\delta_{ij}) \newline
&+ (A_{nmj}\epsilon_{nij} + B_{nmij}\delta_{ij} + )^2 + \delta_{ij} + O(\epsilon_{nij}^3)\newline

&= (2z_{mj}A_{nmj})\epsilon_{nij} + (2z_{mj}A'_{nmj} + A_{nmj}^2)\epsilon_{nij}^2 \newline &+ (2z_{mj}B_{nmj}+1)\delta_{ij} + (2z_{mj}B'_{nmj} + B_{nmj}^2)\delta_{ij}^2 \newline &+ (2z_{mj}C_{nmj}+2A_{nmj}B_{nmj})\epsilon_{nij}\delta_{ij} + O(\epsilon_{nij}^3)\end{align}$$

We can then equate terms in (35) and (41). $\square$

## Precision

I need to talk about the precision of Definition 5, as well as the corollary to take only the linear terms.

## Orbit shifting

Suppose we have calculated an orbit $ik$ (relative to the reference orbit $k$), and we want to calculate an orbit $ji$ relative to the secondary orbit $i$.

If the point $j$ is much closer to $i$ than to $k$, then $|\delta_{ij}| \ll |\delta_{jk}|$, and we would ideally use orbit $i$ to jump forward instead of orbit $j$, since this will give us much more precision.

**Theorem 4:** (Orbit shifting)

For a calculated orbit $ik$ relative to a reference orbit $k$, we can calculate an orbit $ji$ relative to $i$ with

$$
\begin{align}
\epsilon_{mjk} = A_{nmi}(\epsilon_{njk} - \epsilon_{nik}) + B_{nmi}\delta_{ji} + \epsilon_{mik} + O(\epsilon_{nji}^2)
\end{align}
$$

Proof: Using Theorem 1,

$$\begin{align}
\epsilon_{mjk} &= \epsilon_{mji} + \epsilon_{mik}\newline
\end{align}$$

From Definition 5, using only the linear terms,

$$\begin{align}
&= A_{nmi}\epsilon_{nji} + B_{nmi}\delta_{ji} + \epsilon_{mik}+ O(\epsilon_{nji}^2)\newline
\end{align}$$

Using Theorem 1,

$$\begin{align}
&= A_{nmi}(\epsilon_{njk} - \epsilon_{nik}) + B_{nmi}\delta_{ji} + \epsilon_{mik} + O(\epsilon_{nji}^2)\newline
\end{align}$$

$\square$

Theorem 4 tells us that we can use *any* orbit to calculate our jump steps, not just the reference orbit.

## Turning this into an algorithm

Firstly, we'll calculate a high precision reference orbit $k$ using Definition 1.

Each time we calculate steps, we also store an entry for each iteration $n$ consisting of:

* $n$ - which iteration we are jumping from (may be implicit)
* $m$ - which iteration we are jumping to (may be implicit) 
* Terms $A_nmi$ ... $C_nmi$
* $\delta_{ik}$
* $\epsilon_{nik}$ 


