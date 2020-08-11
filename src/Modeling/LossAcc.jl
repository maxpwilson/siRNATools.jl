function generate_accuracy(m::Chain{T}, type::String) where {T}
    if type == "Pearson"
        pearson_coef(m)
    elseif type == "Spearman"
        spearman_coef(m)
    elseif type == "Kendall"
        kendall_coef(m)
    end
end

function pearson_coef(m::Chain{T}) where {T}
    accuracy(x, y) = cor([i[1] for i in m.(x)], y)
end
function spearman_coef(m::Chain{T}) where {T}
    accuracy(x, y) = corspearman([i[1] for i in m.(x)], y)
end
function kendall_coef(m::Chain{T}) where {T}
    accuracy(x, y) = corkendall([i[1] for i in m.(x)], y)
end

function generate_penalty()
    no_penalty()
end
function generate_penalty(m::Chain{T}, alpha::Float64, type::String) where {T}
    if type == "Ridge"
        ridge_penalty(m, alpha)
    else
        no_penalty()
    end
end
function ridge_penalty(m::Chain{T}, alpha::Float64) where {T}
    penalty() = alpha * sum(norm, params(m))
end
function no_penalty()
    penalty() = 0
end
