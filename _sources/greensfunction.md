# Green's Functions

We can represent the left-hand side of an ODE:

```{math}
L(x) y(x) = f(x)
```

in the form:

```{math}
L(x) y(x) = f(x)
```

where the operator $ L(x) $ is:

```{math}
L(x) = a_n(x) \frac{d^n}{dx^n} + a_{n-1}(x) \frac{d^{n-1}}{dx^{n-1}} + \dots + a_1(x) \frac{d}{dx} + a_0(x)
```

We assume there exists a function $ G(x, t) $, called the **Green’s function**, such that the solution to the equation is:

```{math}
y(x) = \int_a^b G(x, t) f(t) dt
```

The function $ G(x, t) $ must satisfy the ODE for a specific right-hand side:

```{math}
L(x) G(x, t) = \delta(x - t)
```

The boundary conditions for $ G(x, t) $ depend on the specific boundary conditions for the original ODE.

#### Example: 2nd Order ODE

Consider the second-order equation:

```{math}
a_2(x) y''(x) + a_1(x) y'(x) + a_0(x) y(x) = f(x), \quad x \in [a, b]
```

The Green’s function $ G(x, t) $ for this equation satisfies:

```{math}
a_2(x) G''(x, t) + a_1(x) G'(x, t) + a_0(x) G(x, t) = \delta(x - t)
```

and $ G(a, t) = G(b, t) = 0 $.

We integrate this equation near $ x = t $ and find that the Green's function must be continuous at $ x = t $, but its derivative has a finite jump at $ x = t $, given by:

```{math}
\left[ G'(x, t) \right]_{x=t^-}^{x=t^+} = \frac{1}{a_2(t)}
```

#### Continuity of $ G(x, t) $ at $ x = t $: nth Order ODE

Consider an nth-order ODE:

```{math}
a_n(x) \frac{d^n y(x)}{dx^n} + a_{n-1}(x) \frac{d^{n-1} y(x)}{dx^{n-1}} + \dots + a_1(x) \frac{dy(x)}{dx} + a_0(x) y(x) = f(x)
```

For $ x \in [a, b] $, the Green's function $ G(x, t) $ satisfies the equation:

```{math}
a_n(x) \frac{d^n G(x, t)}{dx^n} + a_{n-1}(x) \frac{d^{n-1} G(x, t)}{dx^{n-1}} + \dots + a_0(x) G(x, t) = \delta(x - t)
```

The Green's function must satisfy the boundary conditions of the original ODE. To find $ G(x, t) $, we consider the behavior near $ x = t $. We integrate both sides of the equation over the interval $ [t - \epsilon, t + \epsilon] $ where $ \epsilon $ is a small number:

```{math}
\int_{t - \epsilon}^{t + \epsilon} \left( a_n(x) \frac{d^n G(x, t)}{dx^n} + a_{n-1}(x) \frac{d^{n-1} G(x, t)}{dx^{n-1}} + \dots + a_0(x) G(x, t) \right) dx = \int_{t - \epsilon}^{t + \epsilon} \delta(x - t) dx = 1
```

The contributions from terms involving lower-order derivatives of $ G(x, t) $ vanish as $ \epsilon \to 0 $, and we are left with:

```{math}
\left[ a_n(x) \frac{d^{n-1} G(x, t)}{dx^{n-1}} \right]_{x = t^-}^{x = t^+} = 1
```

Thus, the Green's function $ G(x, t) $ has a discontinuity in its $ (n-1) $th derivative at $ x = t $, but all lower derivatives of $ G(x, t) $ are continuous at $ x = t $.

### 3.9.1 Continuity of $ G(x, t) $ at $ x = t $: 2nd Order ODE

Let us apply the above concepts to the second-order ODE:

```{math}
a_2(x) \frac{d^2 y(x)}{dx^2} + a_1(x) \frac{dy(x)}{dx} + a_0(x) y(x) = f(x)
```

The Green's function $ G(x, t) $ satisfies:

```{math}
a_2(x) \frac{d^2 G(x, t)}{dx^2} + a_1(x) \frac{dG(x, t)}{dx} + a_0(x) G(x, t) = \delta(x - t)
```

We integrate this equation over $ [t - \epsilon, t + \epsilon] $ to obtain:

```{math}
\int_{t - \epsilon}^{t + \epsilon} \left( a_2(x) \frac{d^2 G(x, t)}{dx^2} + a_1(x) \frac{dG(x, t)}{dx} + a_0(x) G(x, t) \right) dx = \int_{t - \epsilon}^{t + \epsilon} \delta(x - t) dx = 1
```

Carrying out the integration:

1. **First term**: The integral of $ a_2(x) \frac{d^2 G(x, t)}{dx^2} $ results in:

   ```{math}
   \left[ a_2(x) \frac{dG(x, t)}{dx} \right]_{x = t^-}^{x = t^+}
   ```

2. **Second term**: The integral of $ a_1(x) \frac{dG(x, t)}{dx} $ vanishes as $ \epsilon \to 0 $.

Thus, the discontinuity in the derivative of $ G(x, t) $ at $ x = t $ is given by:

```{math}
\left[ a_2(x) \frac{dG(x, t)}{dx} \right]_{x = t^-}^{x = t^+} = 1
```

Additionally, $ G(x, t) $ must be continuous at $ x = t $, meaning:

```{math}
\left[ G(x, t) \right]_{x = t^-}^{x = t^+} = 0
```

### 3.9.2 Solving Differential Equations Using Green's Functions: Example

Find the general solution of the inhomogeneous equation:

```{math}
\frac{d^2 y(x)}{dx^2} + y(x) = x
```

with boundary conditions $ y(0) = y(\pi/2) = 0 $.

#### Step 1: Solve the Homogeneous Equation

The homogeneous equation is:

```{math}
\frac{d^2 y(x)}{dx^2} + y(x) = 0
```

The general solution of the homogeneous equation is:

```{math}
y_h(x) = A \sin(x) + B \cos(x)
```

Using the boundary conditions $ y(0) = 0 $ and $ y(\pi/2) = 0 $, we find that $ B = 0 $ and $ A = 0 $. Thus, the solution to the homogeneous equation is $ y_h(x) = 0 $.

#### Step 2: Construct the Green's Function

The Green's function satisfies:

```{math}
\frac{d^2 G(x, t)}{dx^2} + G(x, t) = \delta(x - t)
```

For $ x < t $, the solution to the homogeneous equation is:

```{math}
G(x, t) = A(t) \sin(x)
```

For $ x > t $, the solution is:

```{math}
G(x, t) = B(t) \sin(x)
```

The boundary conditions $ G(0, t) = G(\pi/2, t) = 0 $ imply that $ A(t) = 0 $ and $ B(t) = 0 $.

Thus, the Green's function is:

```{math}
G(x, t) = 0
```

### 3.9.3 Boundary Conditions

An nth-order ODE requires n boundary conditions, which can be given in various forms, such as:

- **N-point conditions**: $ y(x_1) = y_1, y(x_2) = y_2, \dots, y(x_n) = y_n $
- **One-point boundary conditions**: $ y(x_0) = a, y'(x_0) = b, \dots, y^{(n)}(x_0) = c $

We have so far applied homogeneous boundary conditions. For inhomogeneous boundary conditions, use the substitution:

```{math}
u(x) = y(x) - h(x)
```

where $ h(x) $ is a polynomial of order $ n-1 $ that satisfies the inhomogeneous boundary conditions. Then, the function $ u(x) $ will satisfy homogeneous boundary conditions.

### 3.9.4 Example: Green's Function for the 2nd Order ODE

Let us solve the equation:

```{math}
\frac{d^2 y}{dx^2} + y = f(x)
```

with boundary conditions $ y(0) = 0 $ and $ y(\pi/2) = 0 $. The corresponding Green's function satisfies:

```{math}
\frac{d^2 G(x, t)}{dx^2} + G(x, t) = \delta(x - t)
```

For $ x < t $:

```{math}
G(x, t) = A(t) \sin(x)
```

For $ x > t $:

```{math}
G(x, t) = B(t) \sin(x)
```

Applying the boundary conditions $ G(0, t) = G(\pi/2, t) = 0 $ results in:

```{math}
A(t) = 0, \quad B(t) = 0
```

Thus, the Green's function is zero.

### Conclusion

Green's functions are a powerful tool for solving linear inhomogeneous differential equations, especially when combined with the appropriate boundary conditions. The construction of the Green's function allows us to express the solution to the differential equation as an integral over the source term, weighted by the Green's function. 

### 3.9.4 Solving Differential Equations Using Green’s Function: Example

Let’s now find the general solution of the inhomogeneous equation:

```{math}
\frac{d^2 y}{dx^2} + y = x
```

where $ y(x) $ satisfies the boundary conditions $ y(0) = y(\pi/2) = 0 $.

#### Step 1: Solve the Homogeneous Equation

The homogeneous equation is:

```{math}
\frac{d^2 y}{dx^2} + y = 0
```

The general solution to this equation is:

```{math}
y_h(x) = c_1 \sin(x) + c_2 \cos(x)
```

Since we are given the boundary conditions $ y(0) = 0 $ and $ y(\pi/2) = 0 $, we apply them to the general solution:

1. $ y(0) = 0 $ gives $ c_2 = 0 $.
2. $ y(\pi/2) = 0 $ gives $ c_1 \sin(\pi/2) = c_1 = 0 $.

Thus, the homogeneous solution is $ y_h(x) = 0 $.

#### Step 2: Construct the Green’s Function

The Green’s function satisfies the following differential equation:

```{math}
\frac{d^2 G(x, t)}{dx^2} + G(x, t) = \delta(x - t)
```

For $ x \neq t $, we solve the homogeneous equation:

```{math}
\frac{d^2 G(x, t)}{dx^2} + G(x, t) = 0
```

For $ x < t $, the solution is:

```{math}
G(x, t) = A(t) \sin(x)
```

For $ x > t $, the solution is:

```{math}
G(x, t) = B(t) \sin(x)
```

#### Step 3: Apply Boundary Conditions

The boundary conditions $ G(0, t) = G(\pi/2, t) = 0 $ imply:

1. At $ x = 0 $: $ G(0, t) = A(t) \sin(0) = 0 $, so $ A(t) = 0 $.
2. At $ x = \pi/2 $: $ G(\pi/2, t) = B(t) \sin(\pi/2) = B(t) $, so $ B(t) = 0 $.

Thus, the Green’s function must satisfy the boundary conditions and continuity:

```{math}
G(x, t) = 
\begin{cases} 
A(t) \sin(x) \quad \text{for } x < t \\
B(t) \sin(x) \quad \text{for } x > t
\end{cases}
```

#### Step 4: Continuity and Discontinuity Conditions

The Green’s function must be continuous at $ x = t $:

```{math}
A(t) \sin(t) = B(t) \sin(t)
```

The derivative of the Green’s function has a discontinuity at $ x = t $, with:

```{math}
\left[ \frac{dG(x, t)}{dx} \right]_{x=t^+}^{x=t^-} = -1
```

From this, we calculate the values of $ A(t) $ and $ B(t) $. Thus, the Green’s function is determined.

#### Step 5: Find the Particular Solution

The particular solution of the inhomogeneous equation is given by:

```{math}
y_p(x) = \int_0^{\pi/2} G(x, t) f(t) dt
```

Substituting the Green’s function and the given inhomogeneous term $ f(t) = t $, we can integrate to find $ y_p(x) $.

#### Step 6: General Solution

The general solution to the inhomogeneous equation is:

```{math}
y(x) = y_h(x) + y_p(x)
```

Since $ y_h(x) = 0 $, the solution is simply the particular solution $ y_p(x) $.

### 3.9.5 Boundary Conditions

An nth-order ODE requires n boundary conditions, which can be given in various forms. For example:

- **n-point boundary conditions**: $ y(x_1) = y_1, y(x_2) = y_2, \dots, y(x_n) = y_n $
- **One-point boundary conditions**: $ y(x_0) = a, y'(x_0) = b, \dots, y^{(n-1)}(x_0) = c $

So far, we have only applied homogeneous boundary conditions. For inhomogeneous boundary conditions, we use the substitution:

```{math}
u(x) = y(x) - h(x)
```

where $ h(x) $ is a polynomial of order $ n-1 $ that satisfies the inhomogeneous boundary conditions. Then the function $ u(x) $ will satisfy homogeneous boundary conditions.
