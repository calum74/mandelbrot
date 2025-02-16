*This is a work in progress, and should be viewed as rough notes*

Next idea: When is it reasonable to drop quadratic terms completely?

# Bilinear approximation

Recall the basic recurrence relation to calculate the distance $\epsilon_i$ from a reference orbit at iteration $i$:

1. $\epsilon_{i+1} = 2z_i\epsilon_i + \epsilon_i^2 + \delta$

If $\epsilon_i$ is "small" (less than $10^{-16}$ or so), then we can drop the $\epsilon_i^2$ term entirely, since it cannot affect the calculation. This gives

2. $\epsilon_{i+1} \approxeq 2z_i\epsilon_i + \delta$

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

Define $\epsilon \in R$ to be the smallest number representable using a floating point number such that if $b<\epsilon a$ then $a \approxeq a+b$.

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

## Translating bilinear terms

Suppose you wanted to calculate a series from a different starting point?

**Theorem:** (Translating bilinear orbits)

For two orbits series $z'$ and $z''$, such that $z'' = z' + \Delta$, where the terms of $z'$ can be estimated using the bilinear equations $\epsilon'_i \approxeq A'_{j,i-j}\epsilon'_j + B'_{j,i-j}\delta'$, and the terms of $z''$ can be estimated using $\epsilon''_i \approxeq A''_{j,i-j}\epsilon''_j + B''_{j,i-j}\delta''$, then

$A''_{j,j-i} = A'_{j,j-i}$

$B''_{j,j-i} = B'_{j,j-i}$

*Proof:*

Consider an orbit $p$ close to $z'$ and $z''$.

$\delta' = p - z'$, $\delta'' = p - z''$

Therefore $\delta'' - \delta' = z''-z' = \Delta$.

At iteration $i$, $p_i \approxeq z'_i + A'_{j,i-j}\epsilon'_j + B'_{j,i-j}\delta' \approxeq z''_i + A''_{j,i-j}\epsilon''j + B''_{j,i-j}\delta''$

Let $\Epsilon_i = \epsilon''_i - \epsilon'_i = z''_i - z'_i$, representing the distance between the orbits $z'$ and $z''$.

$\Epsilon_i = A'_i\Epsilon_j+B'_i\Delta$

1. $\epsilon''_i = A''_i\epsilon''_j + B''_i\delta''$
2. $\epsilon''_i = \Epsilon_i - \epsilon'_i$
3. $=A'_i\Epsilon_j+B'\Delta -A'_i\epsilon'_j - B'_i\delta'$
4. $=A'_i\Epsilon_j+B'\Delta -A'_i(\Epsilon_j - \epsilon''_j) - B'_i(\Delta-\delta'')$
5. $=A'_i\epsilon''_j + B'_i\delta''$
$\square$

This is a bit unexpected, but it means that in a given region, the bilinear coefficients are the same. Unfortunately you don't get anything for free because the translation is only valid over the scope of the original linearization.

## Higher order approximation

What if we instead created a quadratic approximation for skip-forward:

$\epsilon_{i+n} = A_{i,n}\epsilon_i + A'_{i,n}\epsilon_i^2 + B_{i,n}\delta + B'_{i,n}\delta^2 + C_{i,n}\epsilon_i\delta + O(\epsilon_i^3)$

Where $A_{i,n}$ is a coefficient for skipping $n$ steps forward from iteration $i$.

$A_{i,0} = 1$, $A'_{i,0} = B_{i,0} = B'_{i,0} = C_{i,0} = 0$

**Theorem 5** (Computing quadratic terms)

$$A_{i,n+1} = 2z_{i+n}A_{i,n}$$

$$A'_{i,n+1} = 2z_{i+n}A'_{i,n} + A_{i,n}^2$$

$$B_{i,n+1} = 2z_{i+n}B_{i,n}+1$$

$$B'_{i,n+1} = 2z_{i+n}B'_{i,n} + B_{i,n}^2$$

$$C_{i,n+1} = 2z_{i+n}C_{i,n}+2A_{i,n}B_{i,n}$$

Proof:

1. $\epsilon_{i+n+1} = A_{i,n+1}\epsilon_i + A'_{i,n+1}\epsilon_i^2 + B_{i,n+1}\delta + B'_{i,n+1}\delta^2 + C_{i,n+1}\epsilon_i\delta + O(\epsilon_i^3)$

from the definition. We also have

2. $\epsilon_{i+n+1} = 2z_{i+n}\epsilon_{i+n} + \epsilon_{i+n}^2 + \delta$

3. $= 2z_{i+n}(A_{i,n}\epsilon_i + A'_{i,n}\epsilon_i^2 + B_{i,n}\delta + B'_{i,n}\delta^2 + C_{i,n}\epsilon_i\delta + O(\epsilon_i^3)) + (A_{i,n}\epsilon_i + A'_{i,n}\epsilon_i^2 + B_{i,n}\delta + B'_{i,n}\delta^2 + C_{i,n}\epsilon_i\delta + O(\epsilon_i^3))^2 + \delta$

4. $= 2z_{i+n}(A_{i,n}\epsilon_i + A'_{i,n}\epsilon_i^2 + B_{i,n}\delta + B'_{i,n}\delta^2 + C_{i,n}\epsilon_i\delta) + (A_{i,n}\epsilon_i + B_{i,n}\delta)^2 + \delta + O(\epsilon_i^3)$ 

5. $= (2z_{i+n}A_{i,n})\epsilon_i + (2z_{i+n}A'_{i,n} + A_{i,n}^2)\epsilon_i^2 + (2z_{i+n}B_{i,n}+1)\delta + (2z_{i+n}B'_{i,n} + B_{i,n}^2)\delta^2 + (2z_{i+n}C_{i,n}+2A_{i,n}B_{i,n})\epsilon_i\delta + O(\epsilon_i^3)$

We can then equate terms in 1. and 5. $\square$

Higher order terms can give a potentially more accurate way to calculate $\epsilon_i$, but they also give a way to quantify the error. This is a more precise way to put bounds on the validity of BLA.

If $|A'_{i,n}\epsilon_i^2+ B'_{i,n}\delta^2 + C_{i,n}\epsilon_i\delta| \le \epsilon|A_{i,n}\epsilon_i + B_{i,n}\delta|$ then we are good.


# Why didn't it work?

When constructing the tree, I was taking the pessimistic view that it would only be possible to jump forward from *any* point in the region. However, some points in the region will have a much smaller $\epsilon_i$, so could potentially be fast-forwarded much further.

So I need to think about how to calculate these additional terms, and when they can be applied.

Another idea is to set the length of the branch to equal 50% of the escape distance (or the bailout, whichever is the lower).

The for root branch, there is no need to compute all terms because all terms have the same $\epsilon = 0$, and there is no need to extend the terms further.

For secondary branches, we'll compute 50% of the remaining iterations?

Each child branch could have a different branching value.

General problem, for a point $\epsilon_i, \delta$, can we find an orbit with a term $A_{i,n},B_{i,n}$ which will safely allow us to skip forward $n$ spaces?

Wait a minute! Can't we just pick a point directly underneath the point being calculated, and calculate a skip-forward with an epsilon always equal to 0???? Surely that can't work. No because we want to calculate the epsilon.

Think of the whole thing as an orbit-cache. If we calculate one orbit, we'll also calculate the terms, so next time we calculate next to the orbit, we'll be able to use the terms.

Why do we get such bad term utilisation in the secondary branches?

Next idea: Each time we calculate an orbit, also calculate jump-forward terms. When we render an adjacent pixel, there's a good chance we can reuse the last orbit calculated.

# A novel algorithm for computing the Mandelbrot set

This algorithm uses bivariate linear approximation (BLA), but uses the last computed point as a new reference orbit for BLA. Since recent points are close (and can be arranged to always be adjacent), it ensures that we always have a "good" reference orbit, leading to much higher "orbit utilisation" when skipping over iterations.

To acheive this, we realise that the BLA terms of nearby orbits are always the same within their radius of validity. So if we have an orbit at one location, then we can simply reuse it for a different location!

If we have terms $\delta'$, $\epsilon'_{i+n}$, $A_{i,n}$, $A'_{i,n}$, $B_{i,n}$, $B'_{i,n}$, $C_{i,n}$, such that

$\epsilon_{i+n} \approxeq A_{i,n}\epsilon_i + A'_{i,n}\epsilon^2_i + B_{i,n}\delta + B'_{i,n}\delta^2 + C_{i,n}\epsilon_i\delta$

$\epsilon'_{i+n} = z_i - \delta'$ calculated

And we have an orbit we are computing $\delta''$, $\epsilon''_i$, relative to $z_i$, then

$\epsilon''_{i+n} \approxeq A_{i,n}\epsilon_i + A'_{i,n}\epsilon^2_i + B_{i,n}\delta + B'_{i,n}\delta^2 + C_{i,n}\epsilon_i\delta + \epsilon'_{i+n}$

$\delta = \delta'' - \delta'$

$\epsilon_i = \epsilon''_i - \epsilon'_i$

Proof:

$\square$

Input:

* Iteration $i$
* $\delta$ from reference
* $\epsilon_i$ from reference

Algorithm: Using the stack, obtain an entry $(i,n)$ that allows us to jump forward n places. 

Theorem: 

If ...

then

$\epsilon_{i+n} \approxeq $

Compute 

We then test that

$A_{i,n} \epsilon_i + B_{i,n} \delta \approxeq A_{i,n} \epsilon_i + B_{i,n} \delta + A'_{i,n} \epsilon_i^2 + B'_{i,n}\delta^2 + C_{i,n}\epsilon_i\delta$


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


## Delete this section

**Theorem 7:** If we have orbits $i$ and $k$, and we have calculated orbit $ik$, we can calculate orbit $ij$ using:

$$
\begin{align}
A_j &= A_k + 2A_k\epsilon_{ik} + C_k\delta_{ik}\newline
A'_j &= -A'_k\newline
B_j &= B_k + 2B_k\delta_{ik} + C_k\epsilon_{ik}\newline
B'_j &= -B'_k\newline
C_j &= -C_k
\end{align}
$$

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
& +A_{nmk}\epsilon_{nik} - A_{nmk}\epsilon_{nij} + A'_{nmk}\epsilon_{nik}^2 -2A'_{nmk}\epsilon_{nik}\epsilon_{nij} + A'_{nmk}\epsilon_{nij}^2 \newline &+ B_{nmk}\delta_{ik} - B_{nmk}\delta_{ij} + B'_{nmk}\delta_{ik}^2 -2B'_{nmk}\delta_{ik}\delta_{ij} + B'_{nmk}\delta_{ij}^2 \newline &+ C_{nmk}\epsilon_{nik}\delta_{ik} - C_{nmk}\epsilon_{nik}\delta_{ij} - C_{nmk}\epsilon_{nij}\delta_{ik} + C_{nmk}\epsilon_{nij}\delta_{ij} + O(\epsilon_{njk}^3)
\newline

\end{align}$$

Gathering like terms, we see that we have some additional terms in $\epsilon_{ik}\epsilon_{ij}$, $\delta_{ik}\delta_{ij}$, $\epsilon_{ik}\delta_{ij}$ and $\epsilon_{ij}\delta_{ik}$, which we need to put into the relevant term for $ij$.

$$\begin{align}

\epsilon_{mik}&= (A_j - A_k - 2A_k\epsilon_{ik}-C_k\delta_{ik})\epsilon_{nij} \newline
&+ (A'_j+A'_k)\epsilon_{nij}^2 \newline
&+ (B_j-B_k-2B_k\delta_{ik}-C_k\epsilon_{ik})\delta_{ij} \newline
&+ (B'_j+B'_k)\delta_{ij}^2 \newline
&+ (C_j+C_k)\epsilon_{nij}\delta_{ij} \newline

&+ (A_k)\epsilon_{nik} \newline
&+ (A'_k)\epsilon_{nik}^2 \newline
&+ (B_k)\delta_{ik} \newline
&+ (B'_k)\delta_{ik}^2 \newline
&+ (C_k)\epsilon_{nik}\delta_{ik} \newline
\end{align}$$

Equating the terms we get

$$
\begin{align}
A_j &= A_k + 2A_k\epsilon_{ik} + C_k\delta_{ik}\newline
A'_j &= -A'_k\newline
B_j &= B_k + 2B_k\delta_{ik} + C_k\epsilon_{ik}\newline
B'_j &= -B'_k\newline
C_j &= -C_k
\end{align}
$$

$\square$

Now, this feels quite fishy.