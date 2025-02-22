# Orbit resetting

Calum Grant, 22nd February 2025

One of the problems with perturbation theory is when the reference orbit escapes, leaving it too short or too imprecise to be used.

Recall the standard equation for computing $\epsilon_n$, the distance from the reference orbit at iteration $n$:

$$\epsilon_{n+1} = 2z_{n}\epsilon_{n} + \epsilon^2_{n} + \delta$$

The problem is that if the reference orbit $z_{n}$ escapes then we have no way to compute $\epsilon_{n+1}$, and we would need to compute a new reference orbit.

A solution, proposed by [Zhuoran](https://fractalforums.org/index.php?topic=4360.0), is to *reset* the orbit which means that we can just go back to the beginning of the reference orbit. This can be formalized as follows:

**Theorem:** (Orbit resetting)

For $z_n, z'_n, c, \delta, \epsilon_n \in\mathbb{C}, n \ge m \ge 0$, if

$$
\begin{align}
z_0 &= 0 \newline
z_{n+1} &= z_n^2 + c \newline
z'_0 &= 0 \newline
z'_{n+1} &= z'^2_{n} + c + \delta \newline
\epsilon_{m} &= z'_{m} \newline
\epsilon_{n+1} &= 2z_{n-m}\epsilon_{n} + \epsilon^2_{n} + \delta \newline
\end{align}
$$

then

$$\begin{align}z'_n = z_{n-m} + \epsilon_{n}\end{align}$$

Proof: By induction on $n$. (7) is our induction hypothesis.

When $n=m$, 

$$\begin{align}z_{n-m} + \epsilon_n &= z_0 + z'_n \newline &= 0 + z'_n \newline &= z'_n\end{align}$$

from (1) and (5). This establishes (7) in the case that $n=m$.

When $n>m$, assume that (7) holds for $n$. 

$$z'_{n+1} = z'^2_{n} + c + \delta$$

from (4)

$$=(z_{n-m}+\epsilon_n)^2 + c + \delta$$

from our induction hypothesis (7)

$$\begin{align}&=z^2_{n-m} + 2z_{n-m}\epsilon_n + \epsilon_n^2 + c + \delta\newline
&=z^2_{n-m} + c + \epsilon_{n+1}\end{align}$$

from (6)

$$=z_{(n+1)-m} + \epsilon_{n+1}$$

from (2). This establishes that if (7) holds for $n$, then (7) holds for $n+1$. 
$\square$

Note we didn't make use of (3).

(6) and (7) allows us to compute $z'_n$ from an arbitrary iteration $m$ in the reference orbit, which means that we will always have a valid reference orbit because we can adjust $m$ to ensure that our reference iteration $n-m$ is always in the valid range.