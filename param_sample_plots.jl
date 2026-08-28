using Plots
using DataFrames
using Arrow

include("parameter_collection.jl")


ex_name = "Ex2"
pi_name = "Π₇"
prefix = "data/b_pi7_norm_30pct"
results_table = DataFrame(Arrow.Table("$(prefix)_sample_results.arrow"))
params_table = DataFrame(Arrow.Table("$(prefix)_param_samples.arrow"))
fitted_params = Pi_7

##########################
### Plotting Functions ###
##########################

function lyp_vs_rss_scatterplot(lyp, rss, title)
    p = scatter(rss, lyp, title = title,
        markerstrokewidth = 0,
        markersize = 3,
        marker = :circle,
        legend = false,
        framestyle = :box,
        xlabel = "Relative RSS",
        ylabel = "Lyapunov Exp.",
    )
    hline!(p, [0], color = :black, label = "")
end


function get_inds_greaterthan(vector, value)
    return findall(>(value), vector)
end

function get_inds_lessthan(vector, value)
    return findall(<(value), vector)
end

function prob_nonchaotic_given_rss_lt_x(rss_vec, lyp_vec, rss_thresh)
    inds_for_rss_lt_thresh = get_inds_lessthan(rss_vec, rss_thresh)
    lyp_vals = lyp_vec[inds_for_rss_lt_thresh]
    pos_lyp_vals = lyp_vals[lyp_vals .< 0]
    # neg_lyp_vals = lyp_vals[lyp_vals .< 0]
    return length(pos_lyp_vals) / length(lyp_vals)
end

function prob_chaotic_given_rss_lt_x(rss_vec, lyp_vec, rss_thresh)
    inds_for_rss_lt_thresh = get_inds_lessthan(rss_vec, rss_thresh)
    lyp_vals = lyp_vec[inds_for_rss_lt_thresh]
    pos_lyp_vals = lyp_vals[lyp_vals .> 0]
    return length(pos_lyp_vals) / length(lyp_vals)
end

function plot_prob_chaotic(rss_vec, lyp_vec, title; delta = 0.01)
    min_rss = minimum(rss_vec)
    max_rss = maximum(rss_vec)
    delta = delta
    eps_rng = min_rss+delta:delta:max_rss
    prob_chaotic_vec = [
        prob_chaotic_given_rss_lt_x(rss_vec, lyp_vec, ep) for ep in eps_rng
    ]
    p = plot(eps_rng, prob_chaotic_vec, title = title,
        linewidth = 2,
        legend = false,
        framestyle = :box,
        xlabel = "ε",
        ylabel = "Pr( λ>0 | RSS < ε )",
        xguidefontsize = 18,
        ylims = (0,1),
    )
    return p
end


function plot_min_viable_param_hist(rss, param_table, rss_max, title)
    bdd_rss_inds = findall(<(rss_max), rss)
    bdd_rss_params = param_table[:, bdd_rss_inds]
    hist_vec = [
        histogram(Vector(bdd_rss_params[i, :]), title = parameter_name_list_unicode[i],
            label = "",
            xrotation = 45
        )
        for i in 1:12
    ]
    for i in 1:12
        vline!(hist_vec[i], [fitted_params[i]], label = "", color = :red, linewidth = 3)
    end
    p = plot(
        hist_vec..., layout = (3,4), size = (800, 600),
        plot_title = title
    )
    return p
end



###########################
### basic scatter plots ###
###########################

p1 = lyp_vs_rss_scatterplot(results_table.lypexp, results_table.rss, "$(ex_name), $(pi_name)")

####################################
### Pr(lyp > 0 | RSS < eps) plot ###
####################################

p2 = plot_prob_chaotic(results_table.rss, results_table.lypexp, "$(ex_name), $(pi_name)")

#############################
### min viable param data ###
#############################

rss_max = 2.0
rss_max_string = string(round(Int, rss_max*100))
p3 = plot_min_viable_param_hist(results_table.rss, params_table, rss_max, "$(ex_name), $(pi_name)")