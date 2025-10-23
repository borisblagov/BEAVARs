using Pkg
Pkg.activate(@__DIR__)
using Documenter, BEAVARs, TimeSeries, LiveServer, DocumenterTools

makedocs(
    sitename="Documentation",
    pages = [
        "Introduction" => "introduction.md",
        "Models" => [
            "Chan2020minn" => "Chan2020minn.md"
            "Chan2020iniw" => "Chan2020iniw.md"
        ],
        "File library" => [
            "Constructors" => "Constructors.md"
            "Initialization" => "init_functions.md"
        ]
    ]
)
# deploydocs(
#     repo = "github.com/borisblagov/BEAVARs.jl.git",
# )