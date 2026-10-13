using CSV, DataFrames, Distributions, LinearAlgebra, Random, Statistics
using PythonPlot

include("model.jl")
include("reduced.jl")

const ROOT   = dirname(@__DIR__)
const DATA   = joinpath(ROOT, "data", "schultz_data.csv")
const FIGDIR = mkpath(joinpath(ROOT, "figures"))

# --- figure style ---------------------------------------------------------------
# Figures are drawn at print size, with one font size throughout.
const FONT  = 12
const WIDTH = 6.5   # full text width, inches

PythonPlot.rc("font", size=FONT)
PythonPlot.rc("axes", labelsize=FONT, titlesize=FONT, linewidth=1.0)
PythonPlot.rc("axes.spines", top=false, right=false)
PythonPlot.rc("xtick", labelsize=FONT)
PythonPlot.rc("ytick", labelsize=FONT)
PythonPlot.rc("legend", fontsize=FONT, frameon=false, handlelength=1.4,
              handletextpad=0.5, labelspacing=0.3, borderaxespad=0.2, columnspacing=1.0)
PythonPlot.rc("lines", linewidth=1.5)
PythonPlot.rc("figure.constrained_layout", use=true)
PythonPlot.rc("pdf", fonttype=42)

const RECOVERED   = "C0"
const ARRESTED    = "C1"
const DRUG_COLORS = ["tab:blue", "tab:orange", "tab:red"]   # low, intermediate, high D
const UG          = 1e-3                                    # ng/ml -> µg/ml

const TIME_LABEL = "Time (\$1/\\lambda_0\$)"
const DRUG_LABEL = "External drug \$D\$ (µg/ml)"

fate_color(df) = df.lambda[end] < 0.5 ? ARRESTED : RECOVERED

panel_label(ax, s) = ax.set_title(s, loc="left", fontweight="bold", fontsize=FONT + 2)

function save(fig, name)
    path = joinpath(FIGDIR, name)
    fig.savefig(path)
    println("wrote ", path)
end
