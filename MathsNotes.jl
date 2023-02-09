module MathsNotes
    using ForwardDiff, LinearAlgebra, Calculus
    export @uselib, d, d², dⁿ, ∂, κ, lHôpital
    macro uselib(lib)
        return :( d(f::Function, x::Number) = $lib.derivative(f, x) )
    end
    d(f::Function) = x -> d(f, x)
    d²(f::Function, x::Number) = d(d(f), x)
    d²(f::Function) = x -> d²(f, x)
    dⁿ(f::Function, n::Integer, x::Number) = n > 0 ? d(n == 1 ? f : x -> dⁿ(f, n-1, x), x) : throw(DomainError(n, "For n-th derivative, n must be positive"))
    dⁿ(f::Function, n::Integer) = x -> dⁿ(f, n, x)
    ∂(f::Function, x::Number; n::Integer = 1) = (v...) -> d(x -> f(v[begin:n-1]..., x, v[n:end]...), x)
    ∂(f::Function; n::Integer = 1) = x -> ∂(f, x; n)
    κ(R::Real) = 1/R
    κ(x::Function, y::Function, t::Real) = (d(x, t)*d²(y, t) - d(y, t)*d²(x, t))/(d(x, t)^2 + d(y, t)^2)^1.5
    κ(s::Vector{Function}, t::Real) = ([d(s[1], t), d(s[2], t), 0]×[d²(s[1], t), d²(s[2], t), 0])[3]/norm([d(s[1], t), d(s[2], t)])^3
    function lHôpital(f::Function, g::Function, a::Number)
        df = d(f, a); dg = d(g, a)
        isnan(f(a)/g(a)) & !isnan(df) & !isnan(dg) ? isnan(df/dg) ? lHôpital(x -> d(f, x), x -> d(g, x), a) : df/dg : throw(ArgumentError("The conditions of l'Hôpital's rule isn't satisfied for f/g at a.\nf(x): $(f(a)), g(x): $(g(a)), df/dx: $df, dg/dx: $dg"))
    end
end
