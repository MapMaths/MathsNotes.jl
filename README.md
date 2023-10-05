# MathsNotes
> Some notes taken when I am studying Mathematical skills, still updating.

## Highlights
These codes will make it more convenient while programming on Maths (Mainly calculus).

## Usage
### Add dependencies

```Julia
julia> ]
(@v1.10) pkg> add Calculus
(@v1.10) pkg> add ForwardDiff
```

### Download and include the codes.

(Actually I don't want to upload these codes to Julia Packages, since they are too simple)

```Julia
julia> include("Desktop/MathsNotes.jl") # Here I put the codes on Desktop
Main.MathsNotes

julia> using .MathsNotes
```

### Select a library (Calculus/ForwardDiff)

Note: The reason I made the macro is that both those two library have their advantages and disadvantages. Therefore, when you find it wrong using one lib, you can change the library you use at any time.

This macro (`@uselib`) will define (or redefine) the function `d`, which all the other functions depend on.

(As reference, I found that the library Calculus support more kinds of functions, while the library ForwardDiff is more accurate)

```Julia
julia> @uselib Calculus
d (generic function with 2 method)
```

### Now use it!

The codes contains functions providing these features:

Real Analysis: derivative, derivative function, second and higher derivative (and their corresponding functions), partial derivative (function), l'Hospital's rule.

Numeric Analysis: Lagrange's Interpolation, ASR (Adaptive Simpson's Rule)

Harmonic Analysis: CFS (Continuous Fourier Series), Whittaker-Shannon's Interpolation

Differential Geometry: curvature

Elementary Number Theory: Rational Root Theorem.

Note: the derivative function calculator will NOT do any symbolic calculations, as it is just a abbreviation for `x -> d(f, x)`

```Julia
julia> @uselib ForwardDiff
d (generic function with 2 methods)

julia> f = d(x -> x)
#1 (generic function with 1 method)

julia> f(8)
1

julia> d(x -> x, Inf)
1.0

julia> d²(sin, π)
-1.2246467991473532e-16

julia> d²(sin, π/2)
-1.0

julia> dⁿ(sin, 4, π/2)
1.0

julia> ∂((x, y) -> x^2 + y^2, 2)(6) # Same as `∂((x, y) -> x^2 + y^2, 2; n=1)(6)`, which `(6)` is the other variable(s).
4

julia> κ(x -> x - sin(x), x -> 1 - cos(x), π) # Curvature
-0.25

julia> lHôpital(x -> sind(180*x), x -> x^2-1, 1)
-1.5707963267948966

julia> @uselib Calculus
d (generic function with 2 methods)

julia> lHôpital(x -> exp(-1/x^2), x -> x^2, 0)
0.0

julia> rrt([3, -31, 102, -111, -1, 30]) # Finding the rational roots of 30x⁵-x⁴-111x³+102x²-31x+3=0, and return the polynomial formed from the irrational roots of the original polynomial additionally.
(Rational[1//2, 1//3, 1//5], Polynomial(-90//1 + 30//1*x + 30//1*x^2))
```

## Need to improve

1. L'Hospital's Rule calculator sometimes gives different answers using different libs.
