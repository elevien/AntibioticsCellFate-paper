# Fig. 4: trajectories of the three-gene model below (A, B) and above (C, D)
# the critical r_y(0), and for a sweep of r_y(0) (E, F).
include(joinpath(@__DIR__, "..", "src", "common.jl"))

p = Params()
#     x       y    z    r_x    r_y  r_z  d    λ
u0 = [0.0001, 0.1, 1.0, 0.001, 0.6, 0.6, 0.0, 1.0]
with_ry(ry) = [u0[1:4]; ry; u0[6:8]]

fig = figure(figsize=(WIDTH, 6.6))
gs  = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1.5])

# --- A-D -----------------------------------------------------------------------
gene_colors = zip(("x", "y", "z"), ("C0", "C3", "C5"))
titles = ("\$r_y(0)<r_y^*\$", "\$r_y(0)>r_y^*\$")
for (j, ry) in enumerate((0.02, 5.0))
    df = simulate(p, with_ry(ry); tmax=10.0)
    top = fig.add_subplot(gs[0, j-1])
    bot = fig.add_subplot(gs[1, j-1], sharex=top)

    top.plot(df.time, df.lambda, "k-", label="\$\\lambda/\\lambda_0\$")
    top.plot(df.time, df.d ./ p.D, "k--", label="\$d/D\$")
    top.axhline(p.Kd / p.D; color="0.5", ls=":", label="\$K_d/D\$")
    top.set_ylim(0, 1.05)
    top.set_yticks([0, 0.5, 1])
    top.set_title(titles[j])
    top.tick_params(labelbottom=false)

    for (g, col) in gene_colors
        bot.plot(df.time, df[!, g], col, label="\$$g\$")
    end
    for (g, col) in gene_colors
        bot.plot(df.time, df[!, "r_"*g], col * "--", label="\$r_$g\$")
    end
    bot.set_xlabel(TIME_LABEL)

    if j == 1
        top.set_ylabel("Growth, drug")
        bot.set_ylabel("Concentration")
    else
        top.legend(loc="center left", bbox_to_anchor=(1.0, 0.5))
        bot.legend(loc="center left", bbox_to_anchor=(1.0, 0.5), ncol=2)
    end
    panel_label(top, "AC"[j:j])
    panel_label(bot, "BD"[j:j])
end

# --- E, F ----------------------------------------------------------------------
runs = [simulate(p, with_ry(ry); tmax=8.0) for ry in 0.0001:0.2:4]

ax3 = fig.add_subplot(gs[2, 0], projection="3d")
ax2 = fig.add_subplot(gs[2, 1])
ax3.view_init(elev=30, azim=20)
ax3.grid(false)
for df in runs
    col = fate_color(df)
    ax3.plot(df.r_y, df.y, df.lambda, "-", color=col, alpha=0.4, lw=1)
    ax3.plot([df.r_y[end]], [df.y[end]], [df.lambda[end]], ".", color=col)
    ax2.plot(df.y, df.lambda, "-", color=col, alpha=0.4, lw=1)
    ax2.plot([df.y[end]], [df.lambda[end]], ".", color=col)
end
ax3.set_xlabel("\$r_y\$")
ax3.set_ylabel("\$y\$")
ax3.set_zlabel("\$\\lambda\$")
ax3.set_box_aspect((1, 1, 0.8), zoom=0.85)
ax2.plot([], [], "-", color=RECOVERED, label="recovered")
ax2.plot([], [], "-", color=ARRESTED, label="arrested")
ax2.legend(loc="center left", bbox_to_anchor=(1.0, 0.5))
ax2.set_xlabel("\$y\$")
ax2.set_ylabel("Growth rate \$\\lambda\$")
ax2.set_yticks([0, 0.5, 1])
panel_label(ax3, "E")
panel_label(ax2, "F")

save(fig, "fig4.pdf")
