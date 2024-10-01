# Independence of Functions

## Linear Independence of Functions

A set of functions $ y_k(x) $ (where $ k=1,2,\dots,n $) is **linearly dependent** if there exist coefficients $ c_i $ (not all zeros) such that:

```{math}
c_1 y_1(x) + c_2 y_2(x) + \dots + c_n y_n(x) = 0
```

Otherwise, these functions are **linearly independent**.

We really need a way to consider the distinct-ness of any solutions we find or if we switch representation of solutions that they form a fundamental set, e.g. $\sin(x),\, \cos(x)$ could be 
considered to be a fundamental set, as could $e^{ix},\, e^{-ix}$, but if we mix the two $\cos(x),\, e^{ix}$ then this is not a fundamental set and we run the risk of having 
duplicate solutions.  

Lets start with a 2nd order system which has two fundamental solutions $y_1,\, y_2$ and has to satisy initial conditions $y(x_0) = y_0,\, y'(x_0) = y_0$:
```{math}
y(x) = Ay_1(x) + By_2(x) \Rightarrow \,&\, y(x_0) = y_0 = Ay_1(x_0) + By_2(x_0)\\
y'(x) = A{y_1}'(x) + B{y_2}'(x) \Rightarrow \,&\, y'(x_0) = {y_0}' = A{y_1}'(x_0) + B{y_2}'(x_0)
```
This means we have two equations to find two unknowns $A,\,B$, we could therefore to solve for these simultaneously, however a different method is to use a 2x2 matrix:
```{math}
\begin{pmatrix} y_0 \\ {y_0}' \end{pmatrix} = \begin{pmatrix} y_1 & y_2 \\ {y_1}' & {y_2}'\end{pmatrix}\,\begin{pmatrix} A \\ B \end{pmatrix}
```
which means to find $A,\,B$ we can use a matrix inverse:
```{math}
\frac{1}{y_1\,{y_2}' - {y_1}'\,y_2}\begin{pmatrix} {y_2}' & -y_2 \\ -{y_1}' & y_1\end{pmatrix}\,\begin{pmatrix} y_0 \\ {y_0}' \end{pmatrix} = \begin{pmatrix} A \\ B \end{pmatrix}
```
Therefore in order for $A,\,B$ to *exist*, the matrix determinant must be non-zero - this is feature here that tells us if we have a funamental set of solutions:
```{math}
W(y_1,\,y_2)(x) = \begin{vmatrix} y_1 & y_2 \\ {y_1}' & {y_2}'\end{vmatrix}
```
We call this determinant the **Wronskian**.

````{admonition} Worked example
:class: seealso

Show that the functions $y_1 = e^{\lambda x}$ and $y_2 = x\,e^{\lambda x}$ form a fundamental set of solutions.
```{math}
W(y_1,\,y_2)(x) = \begin{vmatrix} y_1 & y_2 \\ {y_1}' & {y_2}'\end{vmatrix} &= 
\begin{vmatrix} e^{\lambda x} & x\,e^{\lambda x} \\ \lambda\,e^{\lambda x} & \lambda\,x\,e^{\lambda x} + e^{\lambda x}\end{vmatrix}\\
&= e^{2\lambda x} (\lambda x + 1 - \lambda x) = e^{2\lambda x} 
```
and clearly for any finite $x$, $W \neq 0$, hence we have a fundamental set.
````

````{admonition} Further worked examples
:class: seealso, dropdown
1\. Show that the functions $y_1 = x^{n_1}$ and $y_2 = x^{n_2}$ form a fundamental set of solutions for $x > 0$

Lets switch to $y_1 = e^{n_1\,\ln(x)},\, y_2 = e^{n_2\,\ln(x)}$

```{math}
W(y_1,\,y_2)(x) = \begin{vmatrix} y_1 & y_2 \\ {y_1}' & {y_2}'\end{vmatrix} &= 
\begin{vmatrix} e^{n_1\,\ln(x)} & e^{n_2\,\ln(x)} \\ \frac{n_1}{x}\,e^{n_1\,\ln(x)} & \frac{n_2}{x}\,e^{n_2\,\ln(x)} \end{vmatrix}\\
&= e^{(n_1+n_2)\,\ln(x)} \left(\frac{n_2}{x} - \frac{n_1}{x}\right) = (n_1 - n_1)\,x^{n_1+ n_2 - 1}
```
and clearly for any finite $x$ and $n_1 \neq n_2$, $W \neq 0$, hence we have a fundamental set.


2\. Show that the functions $y_1 = x^{\lambda + \mu i} = e^{(\lambda + \mu i)\ln (x)}$ and 
$y_2 = x^{\lambda - \mu i}= e^{(\lambda - \mu i)\ln (x)}$ form a fundamental set of solutions, where $x > 0$ 

```{math}
W(y_1,\,y_2)(x) = \begin{vmatrix} y_1 & y_2 \\ {y_1}' & {y_2}'\end{vmatrix} &= 
\begin{vmatrix} e^{(\lambda + \mu i)\ln (x)} & e^{(\lambda - \mu i)\ln (x)} \\ \frac{(\lambda + \mu i)}{x}e^{(\lambda + \mu i)\ln (x)} & \frac{(\lambda - \mu i)}{x}e^{(\lambda - \mu i)\ln (x)}\end{vmatrix}\\
&= e^{(\lambda^2 + \mu^2)\ln(x)}\left(\frac{\lambda - \mu i}{x} - \frac{\lambda + \mu i}{x} \right) \\
&= -\frac{2\mu\,i\,e^{(\lambda^2 + \mu^2)\ln(x)}}{x} = -2\mu\,i\,x^{\lambda^2 + \mu^2-1}
```
and clearly for any finite $x > 0$, $W \neq 0$, hence we have a fundamental set.


````

To determine linear independence, differentiate the above equation $ n-1 $ times to obtain a system of $ n $ linear equations:

```{math}
\begin{aligned}
c_1 y_1(x) + c_2 y_2(x) + \dots + c_n y_n(x) &= 0 \\
c_1 y'_1(x) + c_2 y'_2(x) + \dots + c_n y'_n(x) &= 0 \\
\vdots \\
c_1 y_1^{(n-1)}(x) + c_2 y_2^{(n-1)}(x) + \dots + c_n y_n^{(n-1)}(x) &= 0
\end{aligned}
```

The determinant of the matrix formed by the coefficients in front of $ c_k $ is called the **Wronskian**:

```{math}
W(y_1, y_2, \dots, y_n)(x) = 
\begin{vmatrix}
y_1(x) & y_2(x) & \dots & y_n(x) \\
y'_1(x) & y'_2(x) & \dots & y'_n(x) \\
\vdots & \vdots & \dots & \vdots \\
y_1^{(n-1)}(x) & y_2^{(n-1)}(x) & \dots & y_n^{(n-1)}(x)
\end{vmatrix}
```

The Wronskian $ W(y_1, y_2, \dots, y_n)(x) $ is a function of $ x $.

From linear algebra, we know that if $ W(y_1, y_2, \dots, y_n)(x) \neq 0 $, then the only solution of the system of equations is $ c_1 = c_2 = \dots = c_n = 0 $, meaning the functions $ y_i(x) $ are linearly independent. Conversely, if $ W(y_1, y_2, \dots, y_n)(x) = 0 $, then the functions are linearly dependent.

### Example: $ n=2 $

Consider $ y_1(x) = f(x) $ and $ y_2(x) = g(x) $.

A) Assume that $ f(x) $ and $ g(x) $ are linearly dependent, meaning that there are non-zero coefficients $ c_1 $ and $ c_2 $ such that:

```{math}
c_1 f(x) + c_2 g(x) = 0
```

In this case, the Wronskian $ W(f, g)(x) $ is zero. Indeed, if we differentiate the equation above:

```{math}
c_1 f'(x) + c_2 g'(x) = 0
```

The Wronskian is:

```{math}
W(f, g)(x) = 
\begin{vmatrix}
f(x) & g(x) \\
f'(x) & g'(x)
\end{vmatrix} = f(x) g'(x) - f'(x) g(x)
```

Substituting for $ f(x) = -\frac{c_2}{c_1} g(x) $ and $ f'(x) = -\frac{c_2}{c_1} g'(x) $:

```{math}
W(f, g)(x) = -\frac{c_2}{c_1} g(x) g'(x) + \frac{c_2}{c_1} g'(x) g(x) = 0
```

Thus, the Wronskian is zero, indicating linear dependence.

B) Conversely, if $ W(f, g) \neq 0 $, then $ f(x) $ and $ g(x) $ are linearly independent. Consider the system of linear equations:

```{math}
\begin{aligned}
c_1 f(x) + c_2 g(x) &= 0 \\
c_1 f'(x) + c_2 g'(x) &= 0
\end{aligned}
```

This can be represented as a matrix equation $ Ac = b $, where $ A $ is the matrix formed by the coefficients $ f(x), g(x), f'(x), g'(x) $, and $ b = 0 $. Since the determinant $ \text{det}(A) = W(f, g) \neq 0 $, the only solution is $ c_1 = c_2 = 0 $, confirming linear independence.

## The general solution of an inhomogeneous equation 

```{math}
a_n(x) y^{(n)} + a_{n-1}(x) y^{(n-1)} + \dots + a_2(x) y'' + a_1(x) y' + a_0(x) y = f(x)
```

has the form:

```{math}
y(x) = y_c(x) + y_p(x)
```

where $ y_p(x) $ is a particular solution of the inhomogeneous equation and $ y_c(x) $, called the complementary function, is a linear combination of $ n $ linearly independent functions:

```{math}
y_c(x) = c_1 y_1(x) + c_2 y_2(x) + \dots + c_n y_n(x)
```

where each $ y_k(x) $ is a solution of the homogeneous equation:

```{math}
a_n(x) y^{(n)} + a_{n-1}(x) y^{(n-1)} + \dots + a_2(x) y'' + a_1(x) y' + a_0(x) y = 0
```

To check if the obtained solutions $ y_k(x) $ are linearly independent, we use the **Wronskian**.