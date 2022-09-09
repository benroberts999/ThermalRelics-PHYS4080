# Hints

This is a reasonably long project, with quite a lot of coding. I don't want you to get too caught up in the coding details (though that is definitely part of the skill I want you to learn!).
If you get stuck on the code, please collaborate with each other, and ask me questions. I'd rather you spend time on doing and understanding the physics.
Below are a few hints that you may find helpful.
Also in this repository is a few simple example code-snippets. You are free to use these snippets in your own code; of course, you'll have to adept them for the current problem. If you do use these, make sure you can explain how they work and what they're doing.

Some documentation links you might find useful are included
scipy: libraries for Python:
<https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.odeint.html>
<https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html>
<https://docs.scipy.org/doc/scipy/reference/tutorial/integrate.html>
<https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.brentq.html>
<https://docs.scipy.org/doc/scipy/reference/interpolate.html>

GSL (GNU scientific libraries, for c and c++)
<https://www.gnu.org/software/gsl/doc/html/ode-initval.html>
<https://www.gnu.org/software/gsl/doc/html/integration.html>

## Question 2

* The required equation is difficult to solve numerically. There are a few tricks
* With the default parameters, the older python function: `scipy.integrate.odeint` seems to work better than the newer one (`scipy.integrate.solve_ivp`)
* Your life will be made much easier if you define functions for everything you require: a function for Y_eq, a function for dY/dx, a function to convert Y to Omega etc.
* To increase numerical precision available, you should re-scale Y(x) by a constant _before_ solveing the ODE - then rescale back to get the result. For example:

$$
\begin{align}
G(x) &= \mu Y(x)\\
\frac{dG(x)}{dx} &= \mu^2 Y_{eq}^2(x) - G^2(x)
\end{align}
$$

The equation to solve for Y might then look something like:

```python
def solve_dYdx(x_list, spin_j, mass, sv_x):
    x0 = x_list[0]
    y0 = Yeq(spin_j, x0)
    g0 = mu * y0
    return scipy.integrate.odeint(dGdx, g0, x_list, args=(spin_j, mass, sv_x), tfirst=True) / mu
```

Setting the initial x too small will make the equation very difficult to solve (you should think about why).
You can play around with different starting points to see which works best.
A typical starting point may be x~10. You will probably need to integrate out to at least x ~ few 100.

## Question 3

You will find everything much easier if you code everything using natural units, and only convert back to "real" units for input/output. The only really tricky one is converting $\langle\sigma v\rangle$ from cm$^3$/s to GeV$^{-2}$ -- this was done in the lecture notes (be careful to go the right way!)

To find the region in the $m_\chi$ vs. $\langle\sigma v\rangle$ plane where $\Omega\,h^2$ matches the Planck CMB analysis, you can use a root-finding method. One is provided in python's scipy library [docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.brentq](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.brentq.html)

e.g.,

```python
scipy.optimize.brentq(func, a, b)
```

will return x such that f(x)=0, where $a\leq x \leq b$.

In our case, we want function to be

$$
\Omega_{\rm model}(m_\chi, \langle\sigma v\rangle) - \Omega_{\rm observed}
$$

We might, for example, want to loop through several masses, and optimise this function for $\langle\sigma v\rangle$ for each mass. Don't forget the take account of the 3-sigma Planck error bars!

## Question 4

To find the thermally-averaged cross-section, $\langle\sigma v\rangle_{\rm eff}$, you will have to integrate (average) over $s$.

$$
\langle\sigma v\rangle_{\rm eff} = \int_{4m^2}^\infty ..... \, v(s)\,{\rm d}s
$$

For numerical reasons, it is often easier to do a change of variables, t=1/s, and instead integrate over t. Note that t=0 is likely not valid - it instead suffices to take this as a very small number [e.g., $1/(20\,m^2$)]

$$
\langle\sigma v\rangle_{\rm eff} = \int_{0}^{1/(4m^2)} ..... \, v(t)\frac{1}{t^2}\,{\rm d}t
$$

Python's scipy library provides several functions for performing integrals. One that is often a good fit is `scipy.integrate.quad` [docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.quad](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.quad.html)

Providing you write code for the integrand `integrand_sveff_t(t)`, code for this may look something like:

```python
def sv_eff(x, m, lambda_h):
    smin = 4.0 * pow(m, 2)
    smax = 20.0 * pow(m, 2)
    integral, err = scipy.integrate.quad(
        integrand_sveff_t, 1.0/smax, 1.0/smin, args=(x, m, lambda_h), epsabs=0, epsrel=1.0e-3)
    return integral
```

Finally, you will need to interpolate (and extrapolate) the width data given in the table. Again, our friend scipy comes to the rescue; we can use `scipy.interpolate.CubicSpline` [docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.CubicSpline](https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.CubicSpline.html)
Be very careful about the units (GeV vs MeV). The code may look something like this:

```python
Gamma_h = scipy.interpolate.CubicSpline(
    Q_data_GeV, Gamma_data_GeV, extrapolate=True)
```

Gamma_h will now be a _function_ (of Q), that can be sent around and used in the rest of your code.
