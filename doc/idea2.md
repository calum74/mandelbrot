# BLA2

## Introduction

The BLA2 algorithm improves on existing bivariate linear approximation methods by computing higher order terms, and by utilising a novel(?) technique called *orbit shifting* to always locate an efficient reference orbit for bilinear approximation.

### Notation quick reference

$\epsilon_{nij}$ is the distance between two orbits $i$ and $j$ at iteration $n$.

$\delta_{ij}$ is the distance between two orbits at iteration 1.

$A_{nmi}$ and $B_{nmi}$ are bilinear terms used to jump from iteration $n$ to iteration $m$ in orbit $i$.

## Basic definitions

Firstly let's establish some basic notation. When computing a Mandelbrot set, we need to compute whether each point $c_i \in \mathbb{C}$ *escapes*, which means defining a sequence of complex numbers $z_{n,i}$ calculated as follows:

**Definition 1:** Let $z_{n,i} \in \mathbb{C}$ define the sequence of points in orbit $i$, defined by

$$\begin{align}z_{0,i} &= 0 \newline
z_{n+1,i} &= z_{n,i}^2 + c_i\end{align}$$

for some point $c_i \in \mathbb{C}$.
$\square$

Since we will be referring to multiple orbits, we label each orbit $i,j,k...$.

**Definition 2:** Let $\epsilon_{n,i,j} \in \mathbb{C} = z_{n,i} - z_{n,j}$ to be the distance between two orbits $i$ and $j$ at iteration $n$. $\square$

**Definition 3:** Let $\delta_{i,j} \in \mathbb{C} = c_i - c_j$ to be the distance of two orbits. $\square$

Not to be confused with tensor notation and the Kronecker delta. This notation is based on [Wikipedia notation](https://en.wikipedia.org/wiki/Plotting_algorithms_for_the_Mandelbrot_set).

**Theorem 1:** (Basic properties of $\delta_{i,j}$ and $\epsilon_{n,i,j}$)

$$\begin{align}\epsilon_{n,i,j} &= - \epsilon_{n,j,i}\newline
\epsilon_{n,i,i} &= 0\newline
\epsilon_{n,i,k} &= \epsilon_{n,i,j} + \epsilon_{n,j,k} \newline
\delta_{i,j} &= \epsilon_{1,i,j}\newline
\delta_{i,i} &= 0\newline
\delta_{i,j} &= -\delta_{j,i}\newline
\delta_{i,k} &= \delta_{i,j} + \delta_{j,k}
\end{align}$$

Proof: These follow straightforwardly from Definition 3 and 4.

$$\begin{align}
\epsilon_{n,i,j} &= z_{n,i} - z_{n,j} \newline
&= - (z_{n,j} - z_{n,i}) \newline 
&= - \epsilon_{n,j,i}\newline
\epsilon_{n,i,i} &= z_{n,i} - z_{n,i} \newline
&= 0 \newline
\epsilon_{n,i,k} &= z_{n,i} - z_{n,k} \newline
&=z_{n,i} - z_{n,j} + z_{n,j} - z_{n,k} \newline
&= (z_{n,i} - z_{n,j}) + (z_{n,j} - z_{n,k}) \newline
&= \epsilon_{n,i,j} + \epsilon_{n,j,k}
\end{align}$$

The $\delta_{i,j}$ case follows similiarly, or from $\delta_{i,j}=\epsilon_{1,i,j}$. $\square$

Next we'll show a way to iteratively calculate $\epsilon_{n,i,j}$. Traditional perturbation theory only considers a single reference orbit, but we'll extend this to an arbitrary pair of orbits.

**Theorem 2:** (Iterative properties of $\epsilon_{n,i,j}$)
$$\begin{align}\epsilon_{n+1,i,j} = 2z_{n,j}\epsilon_{n,i,j}+ \epsilon_{n,i,j}^2 + \delta_{i,j
}\newline\end{align}$$
Proof:

$$\begin{align}\epsilon_{n+1,i,j} &= z_{n+1,i} - z_{n+1,j} \newline
 &= (z_{n,i}^2 + c_i) - (z_{n,j}^2 + c_j) \newline
 &= ((z_{n,i} + \epsilon_{n,i,j})^2 + c_i) - (z_{n,j}^2 + c_j) \newline
 &= (z_{n,i}^2 + 2z_{n,i}\epsilon_{n,i,j}+ \epsilon_{n,i,j}^2) - z_{n,j}^2 + (c_i - c_j)\newline
 &= 2z_{n,j}\epsilon_{n,i,j} + \epsilon_{n,i,j}^2 + \delta_{i,j}\end{align}$$

$\square$

Making $i$ and $j$ implicit we get the familiar

$$\begin{align}\epsilon_{n+1} = 2z_{n}\epsilon_{n}+ \epsilon_{n}^2 + \delta\end{align}$$ 

## Polynomial approximation

From now on I'll drop the commas and write $\epsilon_{n,i,j}$ as $\epsilon_{nij}$ except where the commas are needed.

**Definition 5:** (Quadratic approximation)

For orbits $i$ and $j$, the distance between the orbits at iteration $m$ is approximated by

$$\epsilon_{mij} = A_{nmj}\epsilon_{nij} + A'_{nmj}\epsilon_{nij}^2 + B_{nmj}\delta_{ij} + B'_{nmj}\delta_{ij}^2 + C_{nmj}\epsilon_{nij}\delta_{ij} + O(\epsilon_{nij}^3)$$

for some constants $A_{nmj}$, $A'_{nmj}$, $B_{nmj}$, $B'_{nmj}$, $C_{nmj}$. $\square$

Note that

$$\begin{align}A_{nnj}=1, A'_{nnj} = B_{nnj} = B'_{nnj} = C_{nnj} = 0\end{align}$$

Since

$$\begin{align}\delta_{ij} = O(\epsilon_{nij})\end{align}$$

follows by Theorem 2, we'll gather up all higher order terms and represent them as $O(\epsilon_{nij}^3)$.

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

## Orbit shifting

Suppose we have calculated an orbit $ik$ (relative to the reference orbit $k$), and we want to calculate an orbit $ji$ relative to the secondary orbit $i$.

If the point $j$ is much closer to $i$ than to $k$, then $|\delta_{ij}| \ll |\delta_{jk}|$ and we should be able to get much better utilization of the orbit.

**Theorem 4:** (Orbit shifting)

For a calculated orbit $ik$ relative to a reference orbit $k$, we can calculate an orbit $ji$ relative to $i$ with

$$
\begin{align}
\epsilon_{mjk} = A_{nmi}(\epsilon_{njk} - \epsilon_{nik}) + B_{nmi}\delta_{ji} + \epsilon_{mik}
\end{align}
$$

Proof: Using Theorem 1,

$$\begin{align}
\epsilon_{mjk} &= \epsilon_{mji} + \epsilon_{mik}\newline
\end{align}$$

From Definition 5,

$$\begin{align}
&= A_{nmi}\epsilon_{nji} + B_{nmi}\delta_{ji} + \epsilon_{mik}\newline

&= A_{nmi}(\epsilon_{njk} - \epsilon_{nik}) + B_{nmi}\delta_{ji} + \epsilon_{mik}\newline
\end{align}$$

$\square$

What Theorem 4 tells us is that we can use *any* orbit to calculate our jump steps, not just the reference orbit.

## Turning this into an algorithm

