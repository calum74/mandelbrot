# Reboot

It still doesn't work so I'll start from the beginning.

We'll manage a set of jump steps $J$, where $J_{i,n,\delta,\epsilon}$ allows us to jump from iteration $i$ to iteration $i+n$
in a linear step.

## Basic definitions

Firstly let's establish some basic notation. When computing a Mandelbrot set, we need to compute whether each point $c \in \mathbb{C}$ *escapes*, which means defining a sequence of complex numbers $z_n$ calculated as follows:

**Definition 1:** Let $z_n \in \mathbb{C}$ to be a point in a *reference orbit* at $c \in \mathbb{C}$, calculated using the standard Mandelbrot set sequence
$$\begin{align}z_0 &= 0 \newline
z_{n+1} &= z_n^2 + c\end{align}$$
$\square$

We use $n$ and $m$ to refer to *iterations* of the orbit.

In perturbation theory, we only need to calculate one reference orbit at $c$, and then we can calulate orbits *relative* to $c$, at a small distance $\delta$ from the orbit. Since we will be referring to multiple orbits, we label each orbit $i,j,k...$, and correspondingly label each $z_{n,i}$ at a distance $\delta_i$ from $c$.


**Definition 2:** Let $z_{n,i} \in \mathbb{C}$ to be the value of the reference orbit at $c+\delta_i$.
$$\begin{align}z_{0,i} &= 0\newline
z_{n+1,i} &= z_{n,i}^2 + c + \delta_i\end{align}$$
$\square$

We can use $0$ to indicate the special orbit at $c$, so $\delta_0 = 0$.

**Definition 3:** Let $\epsilon_{n,i,j} \in \mathbb{C} = z_{n,i} - z_{n,j}$ to be the distance between two orbits $i$ and $j$ at iteration $n$. $\square$

**Definition 4:** Let $\delta_{i,j} \in \mathbb{C} = \delta_i - \delta_j$ to be the distance of two deltas. $\square$

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

Next we'll show a way to iteratively calculate $\epsilon_{n,i,j}$. Traditional perturbation theory does not use a second $j$ orbit and assumes $\delta_j=0$.

**Theorem 2:** (Iterative properties of $\epsilon_{n,i,j}$)
$$\begin{align}\epsilon_{n+1,i,j} = 2z_{n,j}\epsilon_{n,i,j}+ \epsilon_{n,i,j}^2 + \delta_{i,j
}\newline\end{align}$$
Proof:

$$\begin{align}\epsilon_{n+1,i,j} &= z_{n+1,i} - z_{n+1,j} \newline
 &= (z_{n,i}^2 + c + \delta_i) - (z_{n,j}^2 + c + \delta_j) \newline
 &= ((z_{n,i} + \epsilon_{n,i,j})^2 + c + \delta_i) - (z_{n,j}^2 + c + \delta_j) \newline
 &= ((z_{n,i}^2 + 2z_{n,i}\epsilon_{n,i,j}+ \epsilon_{n,i,j}^2) + c + \delta_i) - (z_{n,j}^2 + c + \delta_j) \newline
 &= 2z_{n,j}\epsilon_{n,i,j} + \epsilon_{n,i,j}^2 + \delta_{i,j}\end{align}$$

$\square$

Making $i$ and $j$ implicit we get the familiar

$$\begin{align}\epsilon_{n+1} = 2z_{n}\epsilon_{n}+ \epsilon_{n}^2 + \delta\end{align}$$ 

## Polynomial approximation

From now on I'll drop the commas and write $\epsilon_{n,i,j}$ as $\epsilon_{nij}$ except where the commas are needed.

Starting with linear approximation, often called *bivariate linear approximation* (BLA), we can create a *linear step* that allows us to jump forward from iteration $n$ to iteration $m$.

**Definition 5:** (Linear approximation)

$$\epsilon_{mij} = A_{nmij}\epsilon_{nij} + B_{nmij}\delta_{ij} + O(\epsilon_{nij}^2)$$

$\square$

Note that $m\ge n$, meaning you can jump forward 0 steps, in which case $A_{nnij}=1, B_{nnij} = 0$. (We can probably jump backwards with this idea as well, but that might not be useful...)

The next theorem allows us to calculate values for the linear terms $A$ and $B$.

**Theorem 3:**

$$\begin{align}
A_{n,m+1,i,j} &= 2z_{mj}A_{nmij}\newline
B_{n,m+1,i,j} &= 2z_{mj}B_{nmij}+1\newline
\end{align}$$

Proof:

$$\begin{align}
\epsilon_{m+1,i,j} &= A_{n,m+1,i,j}\epsilon_{nij} + B_{n,m+1,i,j}\delta_{ij} + O(\epsilon^2)\newline

\epsilon_{m+1,i,j}  &= 2z_{mj}\epsilon_{mij}+ \epsilon_{mij}^2 + \delta_{ij}\newline

&= 2z_{mj}(A_{nmij}\epsilon_{nij} + B_{nmij}\delta_{ij} + O(\epsilon^2))+ O(\epsilon^2) + \delta_{ij}\newline

&= (2z_{mj}A_{nmij})\epsilon_{nij} + (2z_{mj}B_{nmij}+1)\delta_{ij} + O(\epsilon^2)\newline
\end{align}$$

Equating terms in (30) and (33) we get

$$\begin{align}
A_{n,m+1,i,j} &= 2z_{mj}A_{nmij}\newline
B_{n,m+1,i,j} &= 2z_{mj}B_{nmij}+1\newline
\end{align}$$

$\square$

Note we haven't said anything about the size of the $\epsilon^2$ term or when this equation is valid.

Next, what happens to the terms $A$ and $B$ in nearby orbits? This is interesting because we want to be able to relocate orbits.

**Theorem 4:** (Translating a linear orbit)

$$\begin{align}
A_{nmik} & = A_{nmij} = A_{nmjk}\newline
B_{nmik} & = B_{nmij} = B_{nmjk}
\end{align}$$

Proof: We'll apply the properties of $\delta$ and $\epsilon$ from Theorem 1.

$$\begin{align}
\epsilon_{mik} &= \epsilon_{mij} + \epsilon_{mjk}\newline

&= A_{nmik}\epsilon_{nik} + B_{nmik}\delta_{ik} + O(\epsilon^2) \newline

&= A_{nmij}\epsilon_{nij} + B_{nmij}\delta_{ij} + A_{nmjk}\epsilon_{njk} + B_{nmjk}\delta_{jk} + O(\epsilon^2)
\newline

&= A_{nmij}(\epsilon_{nik} - \epsilon_{njk}) + B_{nmij}(\delta_{ik} - \delta_{jk}) + A_{nmjk}\epsilon_{njk} + B_{nmjk}\delta_{jk} + O(\epsilon^2)
\newline

&= A_{nmij}\epsilon_{nik} - A_{nmij}\epsilon_{njk} + B_{nmij}\delta_{ik} - B_{nmij}\delta_{jk} + A_{nmjk}\epsilon_{njk} + B_{nmjk}\delta_{jk} + O(\epsilon^2)
\newline

&= A_{nmij}\epsilon_{nik} + (A_{nmjk} - A_{nmij})\epsilon_{njk} + B_{nmij}\delta_{ik} + (B_{nmjk}- B_{nmij})\delta_{jk}  + O(\epsilon^2)
\newline

\end{align}
$$

Equating terms in (38) and (43), we see that

$$\begin{align}
A_{nmik} & = A_{nmij} = A_{nmjk}\newline
B_{nmik} & = B_{nmij} = B_{nmjk}
\end{align}$$

$\square$

This does not mean that the terms for $A$ and $B$ are constant, but that they are the same within a "local" region (to be determined).

## Higher order terms

Having warmed up on linear approximation, let's add some more terms. This is important because it will allow us to better calculate the "region of validity" for the linear case, where we can look at the higher order terms to see whether the equation is valid.

**Definition 6:** (Quadratic approximation)

For orbits $i$ and $j$, the distance between the orbits at iteration $m$ is approximated by

$$\epsilon_{mij} = A_{nmij}\epsilon_{nij} + A'_{nmij}\epsilon_{nij}^2 + B_{nmij}\delta_{ij} + B'_{nmij}\delta_{ij}^2 + C_{nmij}\epsilon_{nij}\delta_{ij} + O(\epsilon_{nij}^3)$$

For some constants $A_{nmij}$, $A'_{nmij}$, $B_{nmij}$, $B'_{nmij}$, $C_{nmij}$. $\square$

Note that

$$\begin{align}A_{nnij}=1, A'_{nnij} = B_{nnij} = B'_{nnij} = C_{nnij} = 0\end{align}$$

Since

$$\begin{align}\delta_{ij} = O(\epsilon_{nij})\end{align}$$

follows by Theorem 2, we'll gather up all higher order terms and represent them as $O(\epsilon_{nij}^3)$.

**Theorem 6:** (Computation of quadratic terms.)

$$\begin{align}
A_{n,m+1,i,j} &= 2z_{mj}A_{nmij}\newline

A'_{n,m+1,i,j} &= 2z_{mj}A'_{nmij} + A_{nmij}^2\newline

B_{n,m+1,i,j} &= 2z_{mj}B_{nmij}+1\newline

B'_{n,m+1,i,j} &= 2z_{mj}B'_{nmij} + B_{nmij}^2\newline

C_{n,m+1,i,j} &= 2z_{mj}C_{nmij}+2A_{nmij}B_{nmij}\end{align}$$

Proof:

$$\begin{align}

\epsilon_{n,m+1,i,j} &= A_{n,m+1,i,j}\epsilon_{nij} + A'_{n,m+1,i,j}\epsilon_{nij}^2 \newline
&+ B_{n,m+1,i,j}\delta_{ij} + B'_{n,m+1,i,j}\delta_{ij}^2 \newline 
&+ C_{n,m+1,i,j}\epsilon_{nij}\delta_{ij} + O(\epsilon_{nij}^3)\newline

\epsilon_{n,m+1,i,j} &= 2z_{mj}\epsilon_{mij} + \epsilon_{mij}^2 + \delta_{ij}\newline

&= 2z_{mj}(A_{nmij}\epsilon_{nij} + A'_{nmij}\epsilon_{nij}^2 + B_{nmij}\delta_{ij} + B'_{nmij}\delta_{ij}^2 + C_{nmij}\epsilon_{nij}\delta_{ij} + O(\epsilon_{nij}^3)) \newline
&+ (A_{nmij}\epsilon_{nij} + A'_{nmij}\epsilon_{nij}^2 + B_{nmij}\delta_{ij} + B'_{nmij}\delta_{ij}^2 + C_{nmij}\epsilon_{nij}\delta_{ij} + O(\epsilon_{nij}^3))^2 + \delta_{ij}\newline
\end{align}
$$
Gather the higher order terms into $O(\epsilon_{nij}^3)$:
$$\begin{align}
&= 2z_{mj}(A_{nmij}\epsilon_{nij} + A'_{nmij}\epsilon_{nij}^2 + B_{nmij}\delta_{ij} + B'_{nmij}\delta_{ij}^2 + C_{nmij}\epsilon_{nij}\delta_{ij}) \newline
&+ (A_{nmij}\epsilon_{nij} + B_{nmij}\delta_{ij} + )^2 + \delta_{ij} + O(\epsilon_{nij}^3)\newline

&= (2z_{mj}A_{nmij})\epsilon_{nij} + (2z_{mj}A'_{nmij} + A_{nmij}^2)\epsilon_{nij}^2 + (2z_{mj}B_{nmij}+1)\delta_{ij} + (2z_{mj}B'_{nmij} + B_{nmij}^2)\delta_{ij}^2 \newline &+ (2z_{mj}C_{nmij}+2A_{nmij}B_{nmij})\epsilon_{nij}\delta_{ij} + O(\epsilon_{nij}^3)\end{align}$$

We can then equate terms in (53) and (61). $\square$

**Corollary 1:**

$$\begin{align}
A_{nmij} = A_{nmkj} = A_{nmj}\newline
A'_{nmij} = A'_{nmkj} = A'_{nmj}\newline
B_{nmij} = B_{nmkj}= B_{nmj}\newline
B'_{nmij} = B'_{nmkj} = B'_{nmj}\newline
C_{nmij} = C_{nmkj} = C_{nmj}\newline
\end{align}
$$

Proof: By construction of Theorem 6, only the base orbit $j$ is used in the construction of each term. $\square$

Corollary 1 means that we only need to specify a single orbit when describing a term.


**Theorem 7:** Triangle rule for polynomial terms.

Assume we know $\epsilon_{nij}$ and we want to calculate the terms for $\epsilon_{nik}$

$$
A_{nmik} = A_{nmij} + ...
$$

or something like that.

Proof:

From Definition 6:

$$\begin{align}
\epsilon_{mik} &= A_{nmk}\epsilon_{nik} + A'_{nmk}\epsilon_{nik}^2 + B_{nmk}\delta_{ik} + B'_{nmk}\delta_{ik}^2 + C_{nmk}\epsilon_{nik}\delta_{ik}+ O(\epsilon_{nik}^3)\newline
\end{align}$$

Using Theorem 1:

$$\begin{align}
\epsilon_{mik} &= \epsilon_{mij} + \epsilon_{mjk}\newline

&= A_{nmj}\epsilon_{nij} + A'_{nmj}\epsilon_{nij}^2 + B_{nmij}\delta_{ij} + B'_{nmj}\delta_{ij}^2 + C_{nmj}\epsilon_{nij}\delta_{ij} \newline

& +A_{nmk}\epsilon_{njk} + A'_{nmj}\epsilon_{njk}^2 + B_{nmk}\delta_{jk} + B'_{nmk}\delta_{jk}^2 + C_{nmk}\epsilon_{njk}\delta_{jk} + O(\epsilon_{njk}^3)
\newline

\end{align}$$

However by Theorem 2, $\epsilon_{njk} = \epsilon_{nji} + \epsilon_{nik} = \epsilon_{nik} - \epsilon_{nij}$, so

$$\begin{align}
\epsilon_{mik}&= A_{nmj}\epsilon_{nij} + A'_{nmj}\epsilon_{nij}^2 + B_{nmj}\delta_{ij} + B'_{nmj}\delta_{ij}^2 + C_{nmj}\epsilon_{nij}\delta_{ij} \newline
& +A_{nmk}(\epsilon_{nik} - \epsilon_{nij}) + A'_{nmk}(\epsilon_{nik} - \epsilon_{nij})^2 \newline &+ B_{nmk}(\delta_{nik} - \delta_{nij}) + B'_{nmk}(\delta_{ik} - \delta_{ij})^2 \newline &+ C_{nmk}(\epsilon_{nik} - \epsilon_{nij})(\delta_{ik} - \delta_{ij}) + O(\epsilon_{njk}^3)
\newline

\end{align}$$

Expanding and distributing,

$$\begin{align}
\epsilon_{mik}&= A_{nmj}\epsilon_{nij} + A'_{nmj}\epsilon_{nij}^2 + B_{nmj}\delta_{ij} + B'_{nmj}\delta_{ij}^2 + C_{nmj}\epsilon_{nij}\delta_{ij} \newline
& +A_{nmk}\epsilon_{nik} - A_{nmk}\epsilon_{nij} + A'_{nmk}\epsilon_{nik}^2 -2A'_{nmk}\epsilon_{nik}\epsilon_{nij} + A'_{nmk}\epsilon_{nij}^2 \newline &+ B_{nmk}\delta_{nik} - B_{nmk}\delta_{nij} + B'_{nmk}\delta_{ik}^2 -2B'_{nmk}\delta_{ik}\delta_{ij} + B'_{nmk}\delta_{ij}^2 \newline &+ C_{nmk}\epsilon_{nik}\delta_{ik} - C_{nmk}\epsilon_{nik}\delta_{ij} - C_{nmk}\epsilon_{nij}\delta_{ik} + C_{nmk}\epsilon_{nij}\delta_{ij} + O(\epsilon_{njk}^3)
\newline

\end{align}$$

Gathering like terms

$$\begin{align}

\epsilon_{mik}&= (A_{nmj} - A_{nmk})\epsilon_{nij} + A'_{nmj}\epsilon_{nij}^2 + (B_{nmj}-B_{nmk})\delta_{ij} + B'_{nmj}\delta_{ij}^2 + C_{nmj}\epsilon_{nij}\delta_{ij} \newline
& +A_{nmk}\epsilon_{nik} + A'_{nmk}(\epsilon_{nik} - \epsilon_{nij})^2 \newline &+ B_{nmk}\delta_{nik} + B'_{nmk}(\delta_{ik} - \delta_{ij})^2 \newline &+ C_{nmk}(\epsilon_{nik} - \epsilon_{nij})(\delta_{ik} - \delta_{ij}) + O(\epsilon_{njk}^3)
\newline

\end{align}$$

Let's equate the terms for $\epsilon_{nij} \implies A_{nmj} = A_{nmk}$

$$
A_{mnik} = ... \newline
A'_{mnik} = A'_{mnij} + 
C_{nmik} = C_{nmjk}\newline
$$


## Relocating linear terms

Suppose we have $i,\delta,\epsilon_{i,\delta,0},\delta',\epsilon_{i,\delta',0}, A_{i,n,\delta,0}, B_{i,n,\delta,0}$, can we use this to calculate $\epsilon_{i+n,\delta',0}$?

Can we calculate $A_{i,n,\delta',0}, B_{i,n,\delta',0}$ wrt $A_{i,n,\delta,0}, B_{i,n,\delta,0}$

Well, we have $\epsilon_{i+n,\delta,0} \approxeq A_{i,n,\delta,0}\epsilon_{i,\delta,0} + B_{i,n,\delta,0}\delta$

Theorem:

$$\epsilon_{i+n,\delta',0} \approxeq A_{i,n,\delta,0}\epsilon_{i,\delta',0} + B_{i,n,\delta,0}\delta'$$

Proof:

$$\begin{align}
\epsilon_{i+n+1,\delta',0} & \approxeq A_{i,n+1,\delta',0}\epsilon_{i,\delta',0} + B_{i,n+1,\delta',0}\delta' \newline
\epsilon_{i+n+1,\delta',0} &= 2z_{i+n,0}\epsilon_{i+n,\delta',0}+ \epsilon_{i+n,\delta',0}^2 + \delta'\newline
&\approxeq 2z_{i+n,0}(A_{i,n,\delta,0}\epsilon_{i,\delta',0} + B_{i,n,\delta,0}\delta') + \delta'\newline
&= \newline
\end{align}$$
