###  Laplace Transform

Finding solutions to differential equations can be made easier if the equation is converted to an algebraic form using the Laplace transform. The Laplace transform is an integral transformation defined as:

```{math}
L[f(x)] = \int_0^\infty f(x) e^{-px} dx = F(p)
```

- By definition, $ L[f] $ does not depend on the behavior of $ f(x) $ for $ x < 0 $.
- Note that function $ F(p) $ might be undefined for some values of $ p $.
- Laplace transform is a linear operation:

```{math}
L[c_1 f(x) + c_2 g(x)] = c_1 L[f(x)] + c_2 L[g(x)]
```

It is straightforward to find the Laplace transform $ F(p) $ of any given function $ f(x) $. To find the function $ g(x) $ given its Laplace transform $ G(p) $, refer to a table of Laplace transforms.

#### Examples

1. For $ f(x) = a $ (constant):

```{math}
F(p) = \int_0^\infty a e^{-px} dx = \frac{a}{p}, \quad \text{for } \Re(p) > 0
```

2. For $ f(x) = e^{-ax} $:

```{math}
F(p) = \int_0^\infty e^{-ax} e^{-px} dx = \frac{1}{a + p}, \quad \text{for } \Re(a + p) > 0
```

3. For $ f(x) = \cos(ax) $:

```{math}
F(p) = \int_0^\infty \cos(ax) e^{-px} dx = \frac{p}{a^2 + p^2}, \quad \text{for } \Re(p) > | \Im(a) |
```

#### Laplace Transform of Derivatives

Consider Laplace transforms of $ f'(x) $ and $ f''(x) $.

1. $ L[f'] $:

```{math}
L[f'] = \int_0^\infty f'(x) e^{-px} dx = \left[ e^{-px} f(x) \right]_0^\infty + p \int_0^\infty f(x) e^{-px} dx = -f(0) + p L[f]
```

2. $ L[f''] $:

```{math}
L[f''] = \int_0^\infty f''(x) e^{-px} dx = p^2 L[f] - p f(0) - f'(0)
```

#### Solving ODE Using the Laplace Transform

Applying the Laplace transform method is convenient if:

- The variable $ x $ is defined in the interval $ [0, \infty) $ (e.g., coordinate of a point in a semi-infinite metal rod or time).
- The values of the function $ f(x) $ and its derivatives are zero for $ x = 0 $.

For example, consider the differential equation:

```{math}
f''(x) + 4 f'(x) + 4 f(x) = x^2 e^{-2x}
```

with initial conditions $ f(0) = 0 $ and $ f'(0) = 0 $. Taking the Laplace transform of the equation gives:

```{math}
p^2 L[f] + 4 p L[f] + 4 L[f] = L[x^2 e^{-2x}]
```

Given the initial conditions, we simplify:

```{math}
(p + 2)^2 L[f] = L[x^2 e^{-2x}]
```

Using a table of Laplace transforms, we find:

```{math}
L[x^2 e^{-2x}] = \frac{2}{(p + 2)^3}
```

Thus, the transformed equation becomes:

```{math}
L[f] = \frac{2}{(p + 2)^5}
```

By referring to a Laplace transform table, we find that the solution of the ODE is:

```{math}
f(x) = \frac{1}{12} x^4 e^{-2x}
```
