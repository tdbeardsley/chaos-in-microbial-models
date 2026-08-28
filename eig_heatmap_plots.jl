using Plots
using Arrow
using DataFrames



###################
### import Data ###
###################

title_name = "Ex1, Π₅"
data_import_prefix = "data/a_pi5_bnd50_n961"
param_sample_table = DataFrame(Arrow.Table("$(data_import_prefix)_param_samples.arrow"))

xy_vals_table = DataFrame(Arrow.Table("$(data_import_prefix)_xy_vals.arrow"))
xy_vec = Vector(xy_vals_table[1, :])

results_table = DataFrame(Arrow.Table("$(data_import_prefix)_sample_results.arrow"))



##########################
### plotting Functions ###
##########################

function plot_scatter_map_twotone(xy_table, z_vec, title_prefix)
    xy_vec = Vector(xy_table[1, :])
    x_vals = first.(xy_vec)
    y_vals = last.(xy_vec)
    p = scatter(x_vals, y_vals, marker_z = z_vec, title = "$(title_prefix), λ > 0 is orange ",
        color=cgrad([theme_palette(:default)[1], theme_palette(:default)[2]], [0.5]; categorical=true),
        clims = (-5, 5),
        marker = :square,
        markerstrokewidth = 0,
        label = "",
        legend = false,
    )
    return p
end

function plot_scatter_map_spectrum(xy_table, z_vec, title_prefix)
    xy_vec = Vector(xy_table[1, :])
    x_vals = first.(xy_vec)
    y_vals = last.(xy_vec)
    p = scatter(x_vals, y_vals, marker_z = z_vec, title = "$(title_prefix)",
        marker = :square,
        markerstrokewidth = 0,
        label = "",
    )
    return p
end

function grid_to_matrix(x, y, z)
    length(x) == length(y) == length(z) ||
        throw(DimensionMismatch("x, y, z must be the same length"))

    xs = sort(unique(x))
    ys = sort(unique(y))

    xidx = Dict(v => i for (i, v) in enumerate(xs))
    yidx = Dict(v => i for (i, v) in enumerate(ys))

    Z = fill(NaN, length(ys), length(xs))
    for (xv, yv, zv) in zip(x, y, z)
        Z[yidx[yv], xidx[xv]] = zv
    end

    count(!isnan, Z) == length(z) ||
        @warn "grid is incomplete or has duplicate (x,y) pairs"
    return xs, ys, Z
end

function plot_heat_map_twotone(xy_table, z_vec, title_prefix)
    xy_vec = Vector(xy_table[1, :])
    x_vals = first.(xy_vec)
    y_vals = last.(xy_vec)
    xs, ys, z = grid_to_matrix(x_vals, y_vals, z_vec)
    p = heatmap(xs, ys, z, title = "$(title_prefix), λ > 0 is orange ",
        color=cgrad([theme_palette(:default)[1], theme_palette(:default)[2]], [0.5]; categorical=true),
        clims = (-5, 5),
        label = "",
        legend = false,
        framestyle = :grid
    )
    scatter!(p, [0], [0], color = :black, markersize = 5, label = "")
    return p
end

function plot_heat_map_spectrum(xy_table, z_vec, title_prefix)
    xy_vec = Vector(xy_table[1, :])
    x_vals = first.(xy_vec)
    y_vals = last.(xy_vec)
    xs, ys, z = grid_to_matrix(x_vals, y_vals, z_vec)
    p = heatmap(xs, ys, z, title = "$(title_prefix)",
        label = "",
        framestyle = :grid
    )
    scatter!(p, [0], [0], color = :black, markersize = 5, label = "")
    return p
end



##################
### make plots ###
##################

p1 = plot_scatter_map_twotone(xy_vals_table, results_table.lypexp, title_name)
p2 = plot_scatter_map_spectrum(xy_vals_table, results_table.lypexp, title_name)



####
### heatmaps if you want em' ###

p3 = plot_heat_map_twotone(xy_vals_table, results_table.lypexp, title_name)
p4 = plot_heat_map_spectrum(xy_vals_table, results_table.lypexp, title_name)





