using Plots
using Arrow
using DataFrames


include("parameter_collection.jl")



####################
### load in data ###
####################


e1_data = DataFrame(Arrow.Table("data/a_pi5_re_eig_vectors.arrow"))
e2_data = DataFrame(Arrow.Table("data/b_pi6_re_eig_vectors.arrow"))
e3_data = DataFrame(Arrow.Table("data/b_pi7_re_eig_vectors.arrow"))
e4_data = DataFrame(Arrow.Table("data/c_pi8_re2_eig_vectors.arrow"))


e_list = [e1_data, e2_data, e3_data, e4_data]

param_list = [Pi_5, Pi_6, Pi_7, Pi_8]


##### checking the eigenvalues
v1_data = DataFrame(Arrow.Table("data/a_pi5_re_eig_values.arrow"))
v2_data = DataFrame(Arrow.Table("data/b_pi6_re_eig_values.arrow"))
v3_data = DataFrame(Arrow.Table("data/b_pi7_re_eig_values.arrow"))
v4_data = DataFrame(Arrow.Table("data/c_pi8_re2_eig_values.arrow"))

p0 = scatter(ones(12), abs.(Vector(v1_data[1, :])), 
    yscale = :log10, 
    yticks = [10^(-15), 10^(-10), 10^(-5), 1, 10^(5)], 
    xlims = [0, 4],
    xticks = ([1,2,3], ["Ex 1", "Ex 2", "Ex 3"]),
    title = "Eigenvalues",
    markershape = :hline,
    markersize = 30,
    markerstrokewidth = 2,
    markercolor = :black,
    label = "",
    xtickfontsize = 18,
    ytickfontsize = 14,
    titlesize = 14
)
scatter!(p0, 2 .* ones(12), abs.(Vector(v3_data[1, :])), yscale = :log10,
    markershape = :hline,
    markersize = 30,
    markerstrokewidth = 2,
    markercolor = :black,
    label = ""
)
scatter!(p0, 3 .* ones(12), abs.(Vector(v4_data[1, :])), yscale = :log10,
    markershape = :hline,
    markersize = 30,
    markerstrokewidth = 2,
    markercolor = :black,
    label= ""
)

##########################
### Plotting Functions ###
##########################

function eig_comp_plot(rel_vals_list, title)
    ### rel_vals_list has form:
    # rel_values_list_1 = [abs.(e_list[i].v_12) ./ param_list[i] for i in eachindex(param_list)]
    ###

    p = scatter(log10.(rel_vals_list[1]), label = "Ex1, Π5",
        xticks=(1:12, parameter_name_list_unicode),
        xtickfontsize = 16,
        ytickfontsize = 12,
        title = title,
        markersize = 5,
        ylabel = "log₁₀ |νᵢ|/pᵢ",
        yguidefontsize = 18
    )
    # scatter!(p, log10.(rel_vals_list[2]), label = "Ex2, Π6", markersize = 5)
    scatter!(p, log10.(rel_vals_list[2]), label = "Ex2, Π7", markersize = 5)
    scatter!(p, log10.(rel_vals_list[3]), label = "Ex3, Π8", markersize = 5)
    return p
end

function scale_eigenvectors(ev, p)
    r = abs.(ev) ./ p
    a = argmax(r)
    b = r[a]
    se = ev ./ b
    se = sign(maximum(se)).* se
    return se
end



#########################
### largest eig graph ###
#########################
inds_for_comp = [1, 3, 4]

rel_values_list_1 = [abs.(scale_eigenvectors(e_list[i].v_12, param_list[i])) ./ param_list[i] for i in inds_for_comp]

p1 = eig_comp_plot(rel_values_list_1, "1st Eig-value, Eig-vector Comparison, Log Hess")



################################
### second largest eig graph ###
################################

rel_values_list_3 = [abs.(scale_eigenvectors(e_list[i].v_11, param_list[i])) ./ param_list[i] for i in inds_for_comp]

p3 = eig_comp_plot(rel_values_list_3, "2nd Eig-value, Eig-vector Comparison, Log Hess")

