using Documenter, BEAVARs

makedocs(
    sitename="Documentation",
    pages = [
        "index.md",
        "How To" => "Manual.md",
        "Files" => [
            "Constructors" => "Constructors.md"
        ]
    ]
)