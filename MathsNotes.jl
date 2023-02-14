module MathsNotes
    using ForwardDiff, LinearAlgebra, Calculus, Plots
    export @uselib, d, d², dⁿ, ∂, κ, lHôpital, rrt, lag, simpson, asr, rint, cfs, mixcfs, wsi
    macro uselib(lib)
        return :( d(f::Function, x::Number) = $lib.derivative(f, x) )
    end

    # Derivatives
    d(f::Function) = x -> d(f, x)
    d²(f::Function, x::Number) = d(d(f), x)
    d²(f::Function) = x -> d²(f, x)
    dⁿ(f::Function, n::Integer, x::Number) = n > 0 ? d(n == 1 ? f : x -> dⁿ(f, n-1, x), x) : throw(DomainError(n, "For n-th derivative, n must be positive"))
    dⁿ(f::Function, n::Integer) = x -> dⁿ(f, n, x)

    # Partial Derivative
    ∂(f::Function, x::Number; n::Integer = 1) = (v...) -> d(x -> f(v[begin:n-1]..., x, v[n:end]...), x)
    ∂(f::Function; n::Integer = 1) = x -> ∂(f, x; n)

    # Curvature
    κ(R::Real) = 1/R
    κ(x::Function, y::Function, t::Real) = (d(x, t)*d²(y, t) - d(y, t)*d²(x, t))/(d(x, t)^2 + d(y, t)^2)^1.5
    κ(s::Vector{Function}, t::Real) = ([d(s[1], t), d(s[2], t), 0]×[d²(s[1], t), d²(s[2], t), 0])[3]/norm([d(s[1], t), d(s[2], t)])^3

    # L'Hospital's Rule
    function lHôpital(f::Function, g::Function, a::Number)
        df = d(f, a); dg = d(g, a)
        isnan(f(a)/g(a)) & !isnan(df) & !isnan(dg) ? isnan(df/dg) ? lHôpital(x -> d(f, x), x -> d(g, x), a) : df/dg : throw(ArgumentError("The conditions of l'Hôpital's rule isn't satisfied for f/g at a.\nf(x): $(f(a)), g(x): $(g(a)), df/dx: $df, dg/dx: $dg"))
    end

    function finddiv(n::Integer)
        a = []
        for i = 1:Int(floor(√abs(n)))
            n%i == 0 && append!(a, i^2 == abs(n) ? [i] : [i, n÷i])
        end
        a
    end
    function con(a::Vector, n::Union{Integer, Rational})
        s = a[end]
        for i = 1:length(a)-1
            s += a[i]*n^(length(a)-i)
        end
        s
    end
    # Rational Root Theorem
    function rrt(a::Vector)
        s = con(a, 0) == 0 ? [0] : []
        db = finddiv(a[1])
        de = finddiv(a[end])
        for p ∈ de, q ∈ db
            con(a, p//q) == 0 && push!(s, q == 1 ? p : p//q)
            con(a, -p//q) == 0 && push!(s, q == 1 ? -p : -p//q)
        end
        s
    end

    # Lagrange interpolation
    lag(p::Vector{Tuple{T, S}} where {T <: Real, S <: Real}) = x::Real -> begin
           ans = 0
           for i ∈ p
               prod = 1
               for j ∈ p
                   i[1] == j[1] || (prod *= (x - j[1])/(i[1] - j[1]))
               end
               ans += i[2]*prod
           end
           ans
    end

    # Simpson's formula
    simpson(f::Function, l::Real, r::Real) = (r-l)*(f(l)+f(r)+4f((l+r)/2))/6
    # Adaptive Simpson's Rule
    function asr(f::Function, l::Real, r::Real, eps::AbstractFloat, ans::Number, step::Integer)
        mid = (l + r)/2
        fl = simpson(f, l, mid)
        fr = simpson(f, mid, r)
        abs(fl + fr - ans) <= 15 * eps && step < 0 ? fl + fr + (fl + fr - ans) / 15 : asr(f, l, mid, eps / 2, fl, step - 1) + asr(f, mid, r, eps / 2, fr, step - 1)
    end
    # Riemann Integral (via ASR)
    rint(f::Function, l::Real, r::Real) = asr(f, l, r, 1e-7, simpson(f, l, r), 12)

    # Continuous Fourier Series
    cfs(f::Function, x0::Real, T::Real, k::Integer) = 1/T*rint(x -> f(x)*exp(-2π*im*k*x/T), x0, x0+T)
    mixcfs(f::Function, x0::Real, T::Real, lo::Integer, hi::Integer) = x -> begin
        ans = 0
        for i = lo:hi
            ans += exp(2π*im*i*x/T)*cfs(f, x0, T, i)
        end
        ans
    end

    # Whittaker-Shannon's Interpolation
    wsi(f, n, int) = t -> begin
        res = 0
        for i = zip(f, n)
            res += i[1]*sinc(t/int - i[2])
        end
        res
    end
end
