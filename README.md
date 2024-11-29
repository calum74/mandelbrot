# High precision perturbation

Original work by Calum Grant.

After a certain point, say ay 10^-320, regular perturbation implementations are no longer sufficient because `double` numbers are no longer able to hold even the deltas. `long double` numbers are available on certain platforms, but not all.

To fix this, we can scale each delta with a factor M, which is a very small number like 2^-640. Then we rewrite our perturbation equations with d = M.D and e = M.E.

This gives

```
(1)     E_(n+1) = 2.z_n.E_n + M.E_n^2 + D
```

The Taylor series for E_n in terms of D around a central orbit is given by

```
(2)     E_n = A_n.D + B_n.D^2 + C_n.D^3 + O(D^4)
```

We'll assume that M and D are both very small, so we can ignore O(D^4).  Substituting (2) into (1),

```
(3)     E_(n+1) = 2.z_n.(A_n.D + B_n.D^2 + C_n.D^3 + O(D^4)) + M.(A_n.D + B_n.D^2 + C_n.D^3 + O(D^4))^2 + D

                = (2.z_n.A_n + 1).D + (2.z_n.B_n + M.A_n^2).D^2 + (2.z_n.C_n + 2.M.A_n.B_n).D^3 + O(D^4)
```

The Taylor series coefficients of E satisfy the iterative relation

```
(4)     A_(n+1) = 2.Z_n.A_n + 1
        B_(n+1) = 2.z_n.B_n + M.A_n^2
        C_(n+1) = 2.z_n.C_n + 2.M.A_n.B_n
        D_(n+1) = 2.z_n.D_n + 2.M.A_n.D_n + 2.M.b_n.C_n
        ...
```

We can use the adapted equations (1) and (4) in our Mandelbrot calculations to compute with higher precision.
