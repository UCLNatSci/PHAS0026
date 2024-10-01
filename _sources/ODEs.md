# Linear Ordinary Differential Equations

Differential equations crop up all the time in Physics, so it worth us spending some time to learn many different techniques for solving them.  

````{admonition} Some examples
:class: seealso, dropdown

**Newtonian Mechanics**

Motion of a point mass $m$ subjected to an external force $F$:

```{math}
\frac{\mathrm{d}^2 r}{dt^2} = a = \frac{1}{m}F
```

**Quantum Mechanics**

Schrodinger equation:

```{math}
i \frac{\partial \psi}{\partial t} = \hat{H} \psi \quad (i = \sqrt{-1})
```

**Electrical Circuits**

Let $R$, $C$, and $L$ be parameters of the circuit in the figure below:

![Electric Circuit](LCR.png)

The charge on the capacitor $q(t)$ is a function of time $t$. To find the dependence of the current $I$ on time, find the total voltage drop along the circuit:

```{math}
V(t) = L \frac{\mathrm{d}I(t)}{dt} + R I(t) + \frac{q(t)}{C}
```

Since $I(t) = \frac{\mathrm{d}q(t)}{dt}$, after differentiating with respect to $t$, we obtain:

```{math}
L \frac{\mathrm{d}^2 I(t)}{dt^2} + R \frac{\mathrm{d}I(t)}{dt} + \frac{I(t)}{C} = \frac{\mathrm{d}V(t)}{dt}
```

In general, coefficients may also depend on a variable:

```{math}
A(x) \frac{\mathrm{d}^2 y(x)}{\,\mathrm{d}x^2} + B(x) \frac{\mathrm{d}y(x)}{\,\mathrm{d}x} + C(x)y(x) = D(x)
```
````

## Definitions

#### Ordinary Differential Equation (ODE)

An ordinary differential equation (ODE) is a relation between a function of one variable $u(x)$, the variable $x$, and the function's first and higher derivatives with respect to the variable $x$:

```{math}
\frac{\mathrm{d}u}{\,\mathrm{d}x}, \frac{\mathrm{d}^2 u}{\,\mathrm{d}x^2}, \ldots
```

For example:

```{math}
\frac{\mathrm{d}^2 y(x)}{\mathrm{d}t^2} = \frac{1}{m}f(x)
```

#### Partial Differential Equation (PDE)

A partial differential equation (PDE) is a relation between a function of several variables $u(x_1, x_2, \ldots, x_N)$, variables $x_1, x_2, \ldots, x_N$, and the function’s first and higher derivatives with respect to these variables:

```{math}
\frac{\partial u}{\partial x_k}, \frac{\partial^2 u}{\partial x_k \partial x_n}, \ldots
```

#### Linear Ordinary Differential Equation

A linear ordinary differential equation is a relation of the form:

```{math}
a_0(x) y(x) + a_1(x) y'(x) + a_2(x) y''(x) + a_3(x) y'''(x) + \ldots = b(x)
```

This definition is valid for any order of the DE. Coefficients $a_k$ (for $k = 0,1, \ldots$) and $b$ may or may not be functions of $x$.

Examples:

- Linear: $y' + 5 = -6y$
- Non-linear: $yy' + 5 = -6y$, $y'^2 + 5 = 6xy$, $y' + 5 = -\sin(6y)$

## Linear First-Order ODE

The purpose of this section is to obtain the general solution of the equation:

```{math}
y'(x) + P(x) y(x) = Q(x)
```

For simplicity, we will write $y$, $P$, and $Q$ instead of $y(x)$, $P(x)$, and $Q(x)$.

#### Homogeneous First-Order ODE

First, consider the case where $Q(x) = 0$:

```{math}
y' + P y = 0
```

or

```{math}
\frac{\mathrm{d}y}{\,\mathrm{d}x} = -P y
```

This equation is separable:

```{math}
\frac{\mathrm{d}y}{y} = -P \,\mathrm{d}x
```

The solution can be obtained by integrating both parts:

```{math}
\ln(y) = -\int P \,\mathrm{d}x + C
```

where $C$ is a constant. Thus, $y$ can be written as:

```{math}
y(x) = e^{-\int P(x) \,\mathrm{d}x + C} = e^C e^{-\int P(x) \,\mathrm{d}x} = A e^{-\int P(x) \,\mathrm{d}x}
```

where $A = e^C$ is a constant. It is convenient to introduce the function $S(x)$:

```{math}
S(x) = \int P(x) \,\mathrm{d}x
```

and

```{math}
\frac{\mathrm{d}S(x)}{\,\mathrm{d}x} = P(x)
```

Then:

```{math}
y(x) = A e^{-S(x)}
```

or:

```{math}
y(x) e^{S(x)} = A
```

By differentiating both sides of this equation:

```{math}
\frac{\mathrm{d}}{\,\mathrm{d}x} \left( y e^{S(x)} \right) = y' e^{S(x)} + y e^{S(x)} \frac{\mathrm{d}S}{\,\mathrm{d}x} = e^{S(x)} (y' + P y) = 0
```

#### Inhomogeneous First-Order ODE

For the inhomogeneous equation:

```{math}
y'(x) + P(x) y(x) = Q(x)
```

We can use this form of the perfect derivative to rewrite this equation as:

```{math}
\frac{\mathrm{d}}{\,\mathrm{d}x} \left( y e^{S(x)} \right) = e^{S(x)} (y' + P y) = e^{S(x)} Q
```

where $e^{S(x)} = e^{\int p \,\mathrm{d}x}$ is known as the **integrating factor**.  

Integration of both sides with respect to $x$ gives:

```{math}
y(x) e^{S(x)} = \int Q(x) e^{S(x)} \,\mathrm{d}x + C
```

Thus, the general solution of the linear first-order DE is:

```{math}
y(x) = e^{-S(x)} \left( \int Q(x) e^{S(x)} \,\mathrm{d}x + C \right)
```

where $S(x) = \int P(x) \,\mathrm{d}x$.

#### Example

Find the solution of:

```{math}
y' + y = e^x
```

In our notation, $P(x) = 1$ and $Q(x) = e^x$. First, consider the homogeneous equation and find $S(x)$:

```{math}
S(x) = \int P(x) \,\mathrm{d}x = x \quad \text{and} \quad y(x) = A e^{-x}
```

For the inhomogeneous equation:

```{math}
y(x) = e^{-x} \left( \int e^x e^x \,\mathrm{d}x + C \right) = e^{-x} \left( \frac{1}{2} e^{2x} + C \right) = \frac{1}{2} e^x + C e^{-x}
```

## Linear Second-Order Homogeneous Equations

Consider homogeneous equations with constant coefficients:

```{math}
a_2 \frac{\mathrm{d}^2 y}{\,\mathrm{d}x^2} + a_1 \frac{\mathrm{d}y}{\,\mathrm{d}x} + a_0 y = 0
```

Let $D$ denote the differential operation $D = \frac{\mathrm{d}}{\,\mathrm{d}x}$. Then, the above ODE can be rewritten as:

```{math}
a_2 D^2 y + a_1 D y + a_0 y = a_2 (D - \lambda_1)(D - \lambda_2)y = 0
```

where $\lambda_1$ and $\lambda_2$ are roots of the characteristic equation:

```{math}
\lambda^2 + \frac{a_1}{a_2} \lambda + \frac{a_0}{a_2} = 0
```

with respect to $\lambda$. The solutions of this equation are:

```{math}
y_1 = c_1 e^{\lambda_1 x} \quad \text{and} \quad y_2 = c_2 e^{\lambda_2 x}
```

### Roots of the Characteristic Equation

To solve a second-order ODE, we first need to find the roots of the corresponding characteristic equation. This equation is quadratic and may have three types of solutions:

1. Two different real roots $ \lambda_1 \neq \lambda_2 $
2. Complex conjugate roots $ \lambda_{1,2} = a \pm ib $
3. Two identical real roots $ \lambda_1 = \lambda_2 = \lambda $


#### Case 1: $ \lambda_1 \neq \lambda_2 $

The general solution is a superposition of the two individual solutions:

```{math}
y(x) = c_1 e^{\lambda_1 x} + c_2 e^{\lambda_2 x}
```


#### Case 2: $ \lambda_{1,2} = a \pm ib $

For this case, the solution is:

```{math}
y(x) = c_1 e^{\lambda_1 x} + c_2 e^{\lambda_2 x}
```

Since $ \lambda_1 $ and $ \lambda_2 $ are related, we can rewrite this as:

```{math}
y(x) = e^{a x} \left( A e^{ibx} + B e^{-ibx} \right)
```

Using Euler's formula $ e^{\pm i \theta} = \cos \theta \pm i \sin \theta $, we can write the solution as:

```{math}
y(x) = e^{a x} \left( C_1 \cos(bx) + C_2 \sin(bx) \right)
```

Thus, the general solution for this case is:

```{math}
y(x) = e^{a x} \left( C \sin(bx + \phi) \right)
```





#### Case 3: $ \lambda_1 = \lambda_2 $

In this case, the differential equation can be written as:

```{math}
a_2 (D - \lambda)(D - \lambda) y = 0
```

One solution is:

```{math}
y_1(x) = A e^{\lambda x}
```

But there will still be another solution to this ODE, to find it consider the second solution having the form:

```{math}
y_2(x) = y_1(x)\, f(x)
```

where $f(x)$ is a function to be found - we can think of $f(x)$ here is a **bridging** function between the different 
solutions $y_1,\,y_2$.  

If we differentiate, we find:
```{math}
y_2'(x) &= y_1'(x)\,f(x) + y_1(x)\,f'(x) \\
y_2''(x) &= y_1''(x)\,f(x) + y_1'(x)\,f'(x) + y_1'(x)\,f'(x) + y_1(x)\,f''(x)   \\
&= y_1''(x)\,f(x) + 2y_1'(x)\,f'(x) + y_1(x)\,f''(x)
```

So if we substitute into the ODE and then group terms by common factors:
```{math}
f(x)\left[ a_2\,y_1''(x) + a_1\,y_1'(x) + a_0\,y_1(x) \right] 
+ f'(x)\left[ 2a_2\,y_1'(x) + a_1\,y_1(x)\right] + f''(x)\,y_1(x) = 0
```

The first bracket will disappear as $y_1(x)$ is a solution to the homogeneous equation.  The second bracket also disappears since in this case, with a repeated root $\lambda = \frac{-a_1}{2a_2}$ and so:

```{math}
\lambda = \frac {y_1'(x)}{y_1(x)} = \frac{-a_1}{2a_2} \Rightarrow 2a_2\,y_1'(x) + a_1\,y_1(x)
```
This leaves $f''(x)=0$ and so integrating up:
```{math}
f''(x) &= 0\\
f'(x) &= C_1 \\
f(x) &= C_1\,x + C_2
```
Which means the full solution is:
```{math}
y(x) = \left( C_1\,x + C_2\right)e^{\lambda\,x}
```

We call this method of using known solutions to find unknown solutions **Reduction of Order** - typically we apply this to find solutions to homogeneous equations (we will encounter further methods to solve inhomogeneous equations later). 

Alternatively we can think of this as a sourced first order problem:

```{math}
u(x) = (D - \lambda) y(x)
```

Then the equation becomes:

```{math}
(D - \lambda) u(x) = 0
```

From this, we get:

```{math}
u(x) = B e^{\lambda x}
```

Substituting this into $ u = (D - \lambda) y $, we have:

```{math}
B e^{\lambda x} = y' - \lambda y
```

This is a first-order ODE, with the general solution:

```{math}
y(x) = (B x + C) e^{\lambda x}
```

So the general solution in this case is:

```{math}
y(x) = (B x + C) e^{\lambda x}
```



## Varying Coefficients - Euler Equations

Up until now we have considered the case of constant coefficients in the second order differential case:
```{math}
a_2\,y''(x) + a_1\,y'(x) + a_0\,y(x) = f(x) 
```

But like in the first order problem case, there is nothing to stop us allowing varying coefficients:
```{math}
a_2(x)\,y''(x) + a_1(x)\,y'(x) + a_0(x)\,y(x) = f(x) 

```

Whilst the space of functions we can choose for $a_2(x),\, a_1(x),\, a_0(x)$ is vast, lets start by trying to solve problems of the form:
```{math}
ax^2\,y''(x) + bx^1\,y'(x) +c x^0 y(x)=0
```

These are known as **Euler Equations**.

We do so by using an ansatz of the form $y = x^n$, which means that:
```{math}
x^{n}\left( an(n-1) + bn + c\right) = 0
```
A quadratic (which we call the characteristic equation) that we can solve:
```{math}
an^2 + (b-a)n + c = 0 \Rightarrow n = \frac{a-b \pm \sqrt{b^2-2ab + a^2 - 4ac}}{2a}
```
which means there are three distinct cases to solve for here (to begin with assume $x > 0$ to avoid any additional 
complications with complex solutions)

1\. *Two distinct, real roots $n = n_1,\, n_2$*

Since the ODE {eq}`EulerEqn` here is linear, we can find the superposition of solutions:
```{math}
y = A\,x^{n_1} + B\,x^{n_2},\, \qquad n_1,\, n_2 \in \mathbb{R}
```

2\. *One repeated, real root $n$*

We find that solving the characteristic equation, the only root is:
```{math}
n = \frac{a-b}{2a}
```
giving a solution $y_1 = Ax^n$.  But there should be two solutions to the ODE {eq}`EulerEqn`, so we can use an ansatz of 
the form:
```{math}
y_2 = x^{n}\,f(x)
```
to find $y_2$, by solving for $f(x)$.  

Using this ansatz, the ODE takes the form:
```{math}
y &= x^n\,f \\
y' &= nx^{n-1}\,f + x^n\,f' \\
y'' &= n(n-1)x^{n-2}\,f + 2nx^{n-1}\,f' + x^{n}f''\\
\Rightarrow ax^2\,y′′+bx\,y′+cy &= ax^{n+1}\,f'' + (2an+b)\,x^{n+1}\,f' + (an(n-1) + bn + c)\,x^n\,f = 0
```
This can be simplified straight away since we are using the value of $n$ which satisfies the characteristic equation 
$an(n-1) + bn + c = 0$, therefore the equation reduces to:
```{math}
ax^{n+1}\,f'' + (2an+b)\,x^{n+1}\,f' = 0 
```
using the fact that $\displaystyle n = \frac{a-b}{2a} \Rightarrow 2an + b = a$ means:
```{math}
ax^{n+1}\left(xf'' + f'\right) = 0
```
which can be solved as a linear 1st order ODE in $f'$:
```{math}
\frac{\mathrm{d}}{\mathrm{d}x}\left( xf'\right) = 0 \Rightarrow xf' = A 
```
which means the form of $f(x)$ is given by:
```{math}
f(x) = A\ln(x) + B
```
and therefore the solution here is:
```{math}
y = x^n\,\left(A \ln(x) + B\right)
```

3\. *Two complex roots $n = \lambda \pm i\mu$*

We can solve the characteristic equation here, with the proviso $(b-a)^2 - 4ac < 0$:
```{math}
n = \frac{(a-b)\pm i\sqrt{4ac-(b-a)^2}}{2a} = \lambda \pm i \mu,\qquad \lambda,\, \mu \in \mathbb{R}
```
If we look at the solutions here, lets start with $n = \lambda + i \mu$:
```{math}
y_1 = Ax^{\lambda + i\mu} = Ae^{\ln(x^{\lambda + i \mu})} = Ae^{(\lambda + i\mu)\ln(x)} = Ax^{\lambda}\left( \cos(\mu
\ln(x)) + i \sin(\mu\ln(x))\right)
```
where we have used the Euler function to convert this expression into one with trigonometric and power law functions.  We 
could do a similar exercise with $y_2$, the only difference here being a $-$ sign in the middle, so the overall solution $y = y_1 + y_2$ can be written as:
```{math}
y = y_1 + y_2 = x^{\lambda}\left(A\cos(\mu\ln(x)) + B \sin(\mu\ln(x))\right)
```
where $A,\, B \in \mathbb{C}$ in general, but with real boundary conditions they will turn out to be real.  

### Other ranges of $x$

So far we have taken $x>0$ to ensure all of our solutions here are real, however if we look again at the Euler equation, 
taking a coordinate transformation $t = -x$ and defining a function $z$ such that:
```{math}
z(t) = y(x) = y(-t)
```
then this means that:
```{math}
\frac{\mathrm{d}}{\mathrm{d}t} z(t) = \dot{z}(t) &= - \frac{\mathrm{d}}{\mathrm{d}x} y(x) = y'(x)\\
\frac{\mathrm{d}^2}{\mathrm{d}t^2} z(t) = \ddot{z}(t)&= + \frac{\mathrm{d}^2}{\mathrm{d}x^2} y(x) = y''(x)
```
and therefore {eq}`EulerEqn` reads as:
```{math}
a(-t)^2\ddot{z} + b(-t)(-\dot{z}) + cz &= 0 \\
\Rightarrow at^2\ddot{z} + bt\dot{z} + cz = 0
```
This shows that for $t>0$ we can still have real solutions to this ODE, there would follow the same form as those found 
for $>0$ - but this range corresponds to $x<0$!

Looking at how this could change the solutions, say for case (1) with real, distinct roots in $n$:
```{math}
z(t) = At^{n_1} + Bt^{n_1} \Rightarrow y(x) = A(-x)^{n_1} + B(-x)^{n_2}, \qquad x < 0
```

Therefore the argument of these solutions is always $|x|$ irrespective of whether $x<0$ or $x > 0$, therfore we modify our 
solutions to be of the form:

1\. $y = A\,x^{n_1} + B\,x^{n_2} \longrightarrow y = A\,|x|^{n_1} + B\,|x|^{n_2}$

2\. $y = x^n\,\left(A \ln(x) + B\right) \Longrightarrow y = |x|^n\,\left(A \ln|x| + B\right)$

3\. $y = x^{\lambda}\left(A\cos(\mu\ln(x)) + B \sin(\mu\ln(x))\right) \longrightarrow y = |x|^{\lambda}\left(A\cos(\mu\ln|
x|) + B \sin(\mu\ln|x|)\right)$

We notice that $x=0$ is not covered by these solutions as it corresponds to $cy=0$.

By a similar principle, we could define a variable of the form $t = x - s$, with predicably similar behaviour, therefore 
we can also shift the ODe to be of the form:
```{math}
a(x-s)^2y'' + b(x-s)y' + cy = 0
```
giving solutions of the form:
1\. $y = A\,x^{n_1} + B\,x^{n_2} \longrightarrow y = A\,|x-s|^{n_1} + B\,|x-s|^{n_2}$

2\. $y = x^n\,\left(A \ln(x) + B\right) \Longrightarrow y = |x-s|^n\,\left(A \ln|x-s| + B\right)$

3\. $y = x^{\lambda}\left(A\cos(\mu\ln(x)) + B \sin(\mu\ln(x))\right) \longrightarrow y = |x-s|^{\lambda}\left(A\cos(\mu
\ln|x-s|) + B \sin(\mu\ln|x-s|)\right)$

however it is also clear here that this solution cannot contain the point $x=s$.

````{admonition} Worked example
:class: seealso, 
Solve the following ODE:
```{math}
x^2\,y'' + xy' - y = 0
```
Using the ansatz $y = x^n$, then:
```{math}
y(x) &= x^n \\
y'(x) &= nx^{n-1}\\
y''(x) &= n(n-1)x^{n-2}
```
which means that the ODE can be written as:
```{math}
x^{n}\Big(n(n-1)+n - 1 \Big) = 0 \Rightarrow n^2 - 1 = 0
```
which means that $n = \pm 1$ and therefore the solution will be of the form:
```{math}
y = Ax^1 + Bx^{-1}
```
````


````{admonition} Further worked examples
:class: seealso, dropdown
1\. Solve the ODE $x^2y'' - 2y = 0, \qquad x >0$

Using the ansatz $y = x^n$:
```{math}
x^n \Big(n(n-1) - 2\Big) = 0 \Rightarrow n^2 - n - 2 = (n-2)(n+1) = 0
```
which has roots of $n = -1,\, 2$ and hence
```{math}
y = Ax^2 + Bx^{-1}
```

2\. Solve the ODE $x^2\,y'' - 3x\,y' + 4y = 0, \qquad x >0$

Using the ansatz $y = x^n$:
```{math}
x^{n}\Big(n(n-1)-3n+4 \Big) &= 0 \Rightarrow n^2-4n+4 = 0 \\
(n-2)^2 &= 0
```
which has a repeated root of $n=2$, so one of the solutions is $y = Ax^2$.  We can then use:
```{math}
y &= x^2 \,f(x) \\
y' &= 2x\,f + x^2\,f' \\
y'' &= 2\,f + 2x\,f' + 2x\,f' + x^2\,f'' = 2\,f + 4x\,f' + x^2\,f''
```
substituting these results in we find:
```{math}
x^2\,y'' - 3x\,y' + 4y = x^4\,f''+ (4x^3-3x^3)\,f' + (2x^2-6x^2+4x^2)\,f = 0 
```
which means we have to solve:
```{math}
x^4\,f'' + x^3\,f' = 0 \Rightarrow f'' + \frac{1}{x}f' = 0
```
Using the integrating factor method, with $\mu = e^{\int \frac{1}{x}\,\mathrm{d}x} = x$ we find:
```{math}
\frac{\mathrm{d}}{\mathrm{d}x}\left(x\,f' \right) = 0 \Rightarrow x\,f' = C 
```
and so the solution is found as:
```{math}
f' = \frac{C}{x} \Rightarrow f(x) = C\ln(x)+D
```
where $C,\, D$ are constants to be found, so the final solution is:
```{math}
y = x^2(A + B\ln(x))
```

3\. Find the solution to the following differential equation on any interval not containing $x = - 6$:
```{math}
3 ( x + 6 )^2 y'' + 25 ( x + 6 ) y' - 16 y = 0
```
Using the ansatz $y = |x+6|^n$ we find:
```{math}
|x+6|^n\left(3n(n-1) + 25n - 16 \right) = 0 \Rightarrow 3n^2 + 22n - 16 = 0
```
which is solved for $(3n - 2)(n + 8) = 0$ nd so $n = \frac{2}{3}, \,-8$, meaning that the general solution is of the form:
```{math}
y = A|x+6|^{2/3} + B|x+6|^{-8}
```
````