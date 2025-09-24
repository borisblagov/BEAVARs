using Documenter, BEAVARs

makedocs(
    sitename="Documentation",
    pages = [
        "index.md",
        "How To" => "Manual.md",
        "Subsection" => [
            "CPZ2024.md"
        ]
    ]
)