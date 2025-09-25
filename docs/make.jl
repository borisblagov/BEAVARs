using Documenter, BEAVARs

makedocs(
    sitename="Documentation",
    pages = [
        "Overview" => "index.md",
        "Introduction" => "introduction.md",
        "Files" => [
            "Constructors" => "Constructors.md"
        ]
    ]
)