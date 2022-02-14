### Nico Schlömer

# On the conversion from OSA-UCS to CIEXYZ

##### Abstract

_This article revisits Kobayasi's and Yosiki's algorithm for conversion of
OSA-UCS into XYZ cooordinates. It corrects some mistakes on the involved
functions and initial guesses and shows that that hundreds of thousands of
coordinates can be converted in less than a second with full accuracy._

## Introduction

In 1974, [MacAdam](https://doi.org/10.1364/josa.64.001691) published the
definition of the OSA-UCS color space that tries to adhere particularly well to
experimentally measured color distances. It combines work that had been going
on since the late 1940s. One aspect of OSA-UCS is that, while the conversion
from CIEXYZ coordinates into OSA-UCS $`Lgj`$ coordinates is straightforward,
the conversion the other way around is not. In fact, there is no conversion
method that works solely in elementary functions. Apparently, this had not been
a design goal of OSA-UCS although is severely limits the usability of OSA-UCS.

In 2002, [Kobayasi and Yosiki](https://doi.org/10.1117/12.464524) presented an
algorithm for conversion from $`Lgj`$ to $`XYZ`$ coordinates that leverages
Newton's method for solving nonlinear equation systems. Unfortunately, the
article remains vague at important points and also contains false assertions
about the nature of the involved functions.

In 2013, [Cao et al.](https://doi.org/10.1364/josaa.30.001508) compared
Kobayasi's and Yosiki's approach with some other, more complex methods based on
artificial neural networks and found the latter to be superior.

In the present note, the author aims to iron out the inaccuracies in Kobayasi's article
and improves the efficiency of the algorithm.

## The forward conversion

The conversion from CIEXYZ-100 coordinates to OSA-UCS $`Lgj`$ coordinates is
defined as follows:

- Compute $`x`$, $`y`$ coordinates via
  ```math
  x = \frac{X}{X + Y + Z},\quad y = \frac{Y}{X + Y + Z}.
  ```
- Compute $`K`$ and $`Y_0`$ as
  ```math
  \begin{split}
    K &= 4.4934 x^2 + 4.3034 y^2 - 4.276 x y - 1.3744 x - 2.5643 y + 1.8103,\\
    Y_0 &= Y K.
  \end{split}
  ```
- Compute $`L'`$ and $`C`$ as

  ```math
  \begin{align}
      \label{eq:lc}
      L' &= 5.9 \left(\sqrt[3]{Y_0} - \frac{2}{3} + 0.042 \sqrt[3]{Y_0 - 30}\right)\\
      \nonumber
      C &= \frac{L'}{5.9 \left(\sqrt[3]{Y_0} - \frac{2}{3}\right)}.
  \end{align}
  ```

  (Note that $`L'`$ is $`L`$ in [the original article](https://doi.org/10.1364/josa.64.001691).)

- Compute RGB as

  ```math
    \begin{bmatrix}
      R\\
      G\\
      B
    \end{bmatrix}
    =
    M
    \begin{bmatrix}
      X\\
      Y\\
      Z
    \end{bmatrix}
    \quad\text{with}\quad
    M=\begin{bmatrix}
      +0.7990 & 0.4194 & -0.1648\\
      -0.4493 & 1.3265 & +0.0927\\
      -0.1149 & 0.3394 & +0.7170
    \end{bmatrix}.
  ```

- Compute $`a`$, $`b`$ as
  ```math
  \begin{bmatrix}
    a\\
    b
  \end{bmatrix}
  =
  A
  \begin{bmatrix}
    \sqrt[3]{R}\\
    \sqrt[3]{G}\\
    \sqrt[3]{B}
  \end{bmatrix}
  \quad\text{with}\quad
  A = \begin{bmatrix}
    -13.7 & +17.7 & -4\\
    1.7 & +8 & -9.7
  \end{bmatrix}.
  ```
- Compute $`L`$, $`g`$, $`j`$ as
  ```math
  L = \frac{L' - 14.3993}{\sqrt{2}},\quad g = Ca,\quad j = Cb.
  ```

## The backward conversion

This section describes the conversion from the $`Lgj`$ to the $`XYZ`$ coordinates.

Given $`L`$, we can first compute

```math
L' = L \sqrt{2} + 14.3993.
```

Equation~\eqref{eq:lc} gives the nonlinear relationship between $`L'`$ and $`Y_0`$ from
which we will retrieve $`Y_0`$. First set $`t\coloneqq \sqrt[3]{Y_0}`$ and consider

```math
  0 = f(t) \coloneqq {\left(\frac{L'}{5.9} + \frac{2}{3} - t\right)}^3 - 0.042^3 (t^3 - 30).
```

$`f`$ is a monotonically decreasing cubic polynomial (see [figure](#figure)).

Hence, it has exactly one root that can be found using the classical Cardano formula:

- Expand $`f(t) = at^3 + bt^2 + ct + d`$ with

  ```math
  \begin{split}
    &u = \frac{L'}{5.9} + \frac{2}{3},\quad v = 0.042^3,\\
    &a = -(v + 1),\quad  b = 3u,\quad  c = -3u^2, \quad d = u^3 + 30v.
  \end{split}
  ```

- Compute the depressed form $`\tilde{f}(x)=a(x^3 + px + q)`$:

  ```math
  p = \frac{3ac - b^2}{3a^2},\quad q = \frac{2b^3 - 9abc + 27a^2d}{27a^3}.
  ```

- Compute the root as
  ```math
  t = -\frac{b}{3a} + \sqrt[3]{
    -\frac{q}{2} + \sqrt{{\left(\frac{q}{2}\right)}^2 + {\left(\frac{p}{3}\right)}^3}
  }
  + \sqrt[3]{
    -\frac{q}{2} - \sqrt{{\left(\frac{q}{2}\right)}^2 + {\left(\frac{p}{3}\right)}^3}
  }.
  ```
  Note that the expression in the square root, $`{\left(q/2\right)}^2 +
  {\left(p/3\right)}^3`$ is always positive since, as argued above, $`f`$ has
  exactly one root.

---
##### Remark

_Kobayasi and Yosiki find the root of $`f`$ using Newton's method. A good
initial guess here is $`t = \frac{L'}{5.9} + \frac{2}{3}`$ since the second
term in $`f(t)`$, containing $`0.042^3`$, is comparatively small. Indeed it
typically only takes around 10 iterations to converge to machine precision._

_Cardano's method finds the root at once at the expense of computing one square
root and two cube roots. This approach is found to be about 15 times faster._

---

##### Figure

<img src="https://raw.githubusercontent.com/nschloe/colorio/assets/f.svg"/>             |  <img src="https://raw.githubusercontent.com/nschloe/colorio/assets/psi.svg"/>
:-------------------------:|:-------------------------:

_Left: Graph of $`f(t)`$ \eqref{eq:f} for $`L'=25`$. Note that the root is not
in the turning point, but close to it. This is because of the small second term
in $`f`$. Right: Graph of the function $`\phi`$ for $`L`$, $`g`$, $`j`$
computed from $`X=12`$, $`Y=67`$, $`Z=20`$. The singularity is at
$`w\approx 0.59652046418`$. Note that the function has three roots only the
largest of which is of interest._

---

From here, one can compute

```math
Y_0 = t^3,\quad
C = \frac{L'}{5.9 \left(t - \frac{2}{3}\right)},\quad
a = \frac{g}{C},\quad
b = \frac{j}{C}.
```

With $`a`$ and $`b`$ at hand, it is now possible via equation~\eqref{eq:ab} to pin down
$`(\sqrt[3]{R}, \sqrt[3]{G}, \sqrt[3]{B})`$ to only one degree of freedom, $`w`$.
The exact value of $`w`$ will be found by Newton iteration. The function $`\phi(w)`$
of which a root needs to be found is defined as follows.

> Append the matrix $`A`$~\eqref{eq:ab} with a row such that the new $`3\times3`$-matrix
> $`\tilde{A}`$ is nonsingular and solve
>
> ```math
> \begin{bmatrix}
>   a\\
>   b\\
>   w
> \end{bmatrix}
> =
> \tilde{A}
> \begin{bmatrix}
>   \sqrt[3]{R}\\
>   \sqrt[3]{G}\\
>   \sqrt[3]{B}
> \end{bmatrix}
> ```
>
> (Kobayasi, for instance, appends $`[1, 0, 0]`$ which corresponds to setting
> $`w=\sqrt[3]{R}`$.) Then compute the tentative $`\tilde{X}`$, $`\tilde{Y}`$, $`\tilde{Z}`$
> via~\eqref{eq:m} and further get the corresponding tentative $`\tilde{Y}_0`$
> from~\eqref{eq:KY0}. Then $`\phi(w) = \tilde{Y}_0(w) - Y_0`$.

If the difference between $`\tilde{Y}_0(w)`$ and $`Y_0`$ from~\eqref{eq:gather}
is 0, the correct $`w`$ has been found. Kobayasi states the function $`\phi`$
is "monotone increasing, convex downward, and smooth". Unfortunately, none of
this is true (see [figure](#figure)). In fact, the function has a
singularity at $`w`$ chosen such that the computed tentative $`\tilde{X}`$,
$`\tilde{Y}`$, $`\tilde{Z}`$ sum up to 0 while the individual values of
$`|\tilde{X}|, |\tilde{Y}|, |\tilde{Z}| > 0`$. This happens if the tentative
$`[R, G, B]`$ is orthogonal on $`[1,1,1] M^{-1}`$.

Fortunately, it seems that the function is indeed convex to the right of the
singularity. Newton's method will hence find the correct (largest) root if the
initial guess $`w_0`$ is chosen larger than the root. Since $`w`$ corresponds to
$`\sqrt[3]{R}`$, it is reasonable to chose $`w_0`$ to be the maximum possible value
that $`\sqrt[3]{R}`$ can take, namely that corresponding to $`X=Y=100`$, $`Z=0`$
(see~\eqref{eq:m}), $`w_0=\sqrt[3]{79.9 + 41.94}\approx 4.9575`$.

---
##### Remark

_[Cao et al.](https://doi.org/10.1364/josaa.30.001508) found that the
conversion to from $`Lgj`$ to $`XYZ`$ takes so long that alternative methods
need to be researched. They even find that the Newton iterations sometimes do
not converge, or find the correct result only to a few digits of accuracy. The
author cannot confirm these observations. The computation of hundreds of
thousands of coordinates at once merely takes a second of computation time on a
recent computer ([figure](#figure-1))._

_To achieve this speed, it is important to vectorize all computation, i.e., not to
perform the conversion for each $`Lgj`$-tuple individually one after another, but to
perform all steps on the array. This also means to perform the Newton iteration on all
tuples until the last one of them has converged successfully, even if some already
converge in the first step. The redundant work inflicted by this approach is far
outweighed by the advantages of vectorization._

_All code is published as open-source in [colorio](https://github.com/nschloe/colorio)._

---

##### Figure

<img src="https://raw.githubusercontent.com/nschloe/colorio/assets/speed-absolute.svg"/>             |  <img src="https://raw.githubusercontent.com/nschloe/colorio/assets/speed-relative.svg"/>
:-------------------------:|:-------------------------:

_Computation speed for arrays of $`Lgj`$ values measured with
[colorio](https://github.com/nschloe/colorio). Left: Comparison with CIELAB and CIECAM02. The
conversion of several hundred thousand $`Lgj`$ values takes about 1 second.
Right: Computation speed relative to the evaluation of the cubic root. For
large arrays, the conversion to $`XYZ`$ is about as costly as the evaluation of
35 cubic roots._

---
