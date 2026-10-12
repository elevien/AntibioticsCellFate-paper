using CSV, DataFrames, Distributions, LinearAlgebra, Random, Statistics
using PythonPlot

include("model.jl")
include("reduced.jl")

const ROOT   = dirname(@__DIR__)
const DATA   = joinpath(ROOT, "data", "schultz_data.csv")
const FIGDIR = mkpath(joinpath(ROOT, "figures"))

PythonPlot.rc("font", size=13)
PythonPlot.rc("axes", labelsize=15, titlesize=15, linewidth=1.1)
PythonPlot.rc("xtick", labelsize=12)
PythonPlot.rc("ytick", labelsize=12)
PythonPlot.rc("legend", fontsize=12)
PythonPlot.rc("axes.spines", top=false, right=false)

function save(fig, name)
    path = joinpath(FIGDIR, name)
    fig.savefig(path, bbox_inches="tight")
    println("wrote ", path)
end

const UG = 1e-3   # ng/ml -> µg/ml

panel_label(ax, s; pad=2.5) =
    ax.text(0.04, 0.96, s, transform=ax.transAxes, ha="left", va="top", fontsize=20,
            bbox=Dict("facecolor" => "white", "edgecolor" => "none", "pad" => pad))
