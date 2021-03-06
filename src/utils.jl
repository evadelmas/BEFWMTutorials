using Plots

"""
This function is used internally by the webplot function, it just transform 
an interaction matrix into a list of interaction.
"""
function matrix_to_list(A)
    idx = findall(A .== 1)
    from_sp = (i->i[1]).(idx)
    to_sp = (i->i[2]).(idx)
    return hcat(from_sp, to_sp)
end

"""
This function is also used internally by the webplot function. 
"""
function invert(v, S)
    1 .+ (S .- v)
end

"""
Plot the interaction matrix :tada: with consumer as either rows or columns :tada: 
whichever you prefer! 

    `plot(A, consasrow)`

- use `plot(A, true)` to have consumers in rows
- use `plot(A, false)` to have consumers in columns
"""
function webplot(A::Array{Int64,2}; consasrow::Bool = true)
    S = size(A,1)
    if !consasrow
        floatA = Float64.(A)
        listA = matrix_to_list(floatA')
        plt = scatter(listA[:,2], invert(listA[:,1], S)
            , ms = 3, c = :black
            , size = (300,300), leg = false
            , framestyle = :box, ticks = ([1.5:1:S;], repeat([""], S)), foreground_color_axis = :white
            , xlims = (0.5,S+0.5), ylims = (0.5,S+0.5))    
        plot!([1,S], [S,1], c = :black, linestyle = :dash)
        xlabel!("consumers")
        ylabel!("resources")
    else
        floatA = Float64.(A)
        listA = matrix_to_list(floatA)
        plt = scatter(listA[:,2], invert(listA[:,1],S)
        , ms = 3, c = :black
        , size = (300,300), leg = false
        , framestyle = :box, ticks = ([1.5:1:S;], repeat([""], S)), foreground_color_axis = :white
        , xlims = (0.5,S+0.5), ylims = (0.5,S+0.5))    
        plot!([1,S], [S,1], c = :black, linestyle = :dash)
        ylabel!("consumers")
        xlabel!("resources")
    end
    return plt
end

"""
Plot the interaction matrix :tada: with consumer as either rows or columns :tada: 
whichever you prefer! 

    `plot(out, consasrow, colorextinct)`

- `out` is the BEFWM simulation output
- `consasrow` (Bool) specifies whether we want consumers as rows (true) or column (false)
- `colorextinct` (Bool) specifies whether we want to shade the extinct species
"""
function webplot(out::Dict{Symbol,Any}; consasrow::Bool = true, colorextinct::Bool = false)
    A = out[:p][:A]
    S = size(A,1)
    id_alive = trues(S)
    id_alive[out[:p][:extinctions]] .= false
    wp = webplot(A; consasrow = consasrow)
    if colorextinct
        ii = invert(i,S)
        for i in out[:p][:extinctions]
            if consasrow
                shp1 = Shape([0.5, S+0.5, S+0.5, 0.5], [ii-0.5, ii-0.5, ii+0.5, ii+0.5])
                shp2 = Shape([i-0.5, i+0.5, i+0.5, i-0.5], [0.5, 0.5, S+0.5, S+0.5])    
            else
                shp1 = Shape([i-0.5, i-0.5, i+0.5, i+0.5], [0.5, S+0.5, S+0.5, 0.5])
                shp2 = Shape([0.5, 0.5, S+0.5, S+0.5], [ii-0.5, ii+0.5, ii+0.5, ii-0.5])    
            end
            plot!(shp1, c = :grey, opacity=.3, lw = 0)
            plot!(shp2, c = :grey, opacity=.3, lw = 0)
        end
    end
    return wp
end

"""
Updates the interaction matrix based on the vector of extinct species 

    `updateA(out)`

- `out` is the BEFWM simulation output

Returns a matrix of dimension S where S is the number of persistent species. This function also
checks for disconnected consumers. If it detect any, it will print a message and return a vector
with the list of disconnected species instead of the updated matrix.
"""
function updateA(out)
    A = out[:p][:A]
    S = size(A,1)
    id_alive = trues(S)
    id_alive[out[:p][:extinctions]] .= false
    Anew = A[id_alive, id_alive]
    #check for status change (consumers can't become producers)
    id_th_prod = sum(A, dims = 2) .== 0
    alive_prod = findall((id_th_prod) .& (id_alive))
    alive_prod = [i[1] for i in alive_prod]
    discosp = []
    for i in alive_prod
        if i ??? findall(out[:p][:is_producer])
            println("/!\\ disconnected consumer $i identified as producer")
            push!(discosp, i)
        end
    end
    toreturn = length(discosp) == 0 ? Anew : discosp
    return toreturn
end