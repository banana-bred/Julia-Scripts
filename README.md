# Julia-Scripts
A small collection of Julia utility scripts and their associated wrappers that I use.

## Convolving data with a Gaussian or Cauchy (Lorentz) distribution
For two functions $f, g : [0, \infty) \to \mathbb{R}, we can define their convolution as
```math
(f*g)(x) \equiv \int\limits_{-\infty}^{\infty} f(x')g(x-x')dx'
```
The assumption here is that the data to convolve represents the output of some function that is zero for non-positive arguments.
The lower integration limit can then be truncated to 0.
```math
(f*g)(x) \equiv \int\limits_{0}^{\infty} f(x')g(x-x')dx'
```
I made this for functions $f$ being scattering cross section functions $\sigma(E)$, where $E$ is some positive energy.
There are two posibilities coded for $g$

### Gaussian Distribution
The gaussian distribution is given by
```math
e^{-(x'-x)^2 / (2\gamma^2)}
```
Cross sections $f(x)$ can be convolved with this distribution:
```math
    \tilde f(x)
    =
    \frac{
        \int\limits_{0}^{\infty} dx' \sigma(x') e^{-(x - x')^2 / (2\gamma^2)}
    }{
        \int\limits_{0}^{\infty} dx' e^{-(x - x')^2 / (2\gamma^2)}
    }
    =
    frac{1}{\gamma\sqrt{2\pi}}
    \int\limits_{0}^{\infty} dx' \sigma(x') e^{-(x - x')^2 / (2\gamma^2)}
```
In implementation, the upper integration limit is taken to be the last available value of $x$.
This should not be an issue if the Gaussian width, $\gamma$, is small enough compared to this.








...
Inline math test $\sigma = \sim 7$

Display test
```math
\alpha \to \sigma
```

## Convolve scattering cross sections with a Maxwell-Boltzmann distribution
...

## Fitting interatomic / intermolecular potentials
...

## Using the scripts
The Julia scripts can be run via the command line or with the help of one of the wrapper scripts.
Each wrapper script looks for its Julia script in the folder given by the environmental variable `JULIABIN`.
