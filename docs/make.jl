using Documenter
using ERGMUserterms

DocMeta.setdocmeta!(ERGMUserterms, :DocTestSetup, :(using ERGMUserterms); recursive=true)

makedocs(
    sitename = "ERGMUserterms.jl",
    modules = [ERGMUserterms],
    authors = "Statistical Network Analysis with Julia",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://Statistical-network-analysis-with-Julia.github.io/ERGMUserterms.jl",
        edit_link = "main",
    ),
    repo = Documenter.Remotes.GitHub("Statistical-network-analysis-with-Julia", "ERGMUserterms.jl"),
    pages = [
        "Home" => "index.md",
        "Getting Started" => "getting_started.md",
        "User Guide" => [
            "Term Interface" => "guide/term_interface.md",
            "Templates and Examples" => "guide/templates.md",
            "Validation and Testing" => "guide/validation.md",
            "Benchmarking" => "guide/benchmarking.md",
        ],
        "API Reference" => [
            "Types" => "api/types.md",
            "Validation" => "api/validation.md",
            "Utilities" => "api/utilities.md",
        ],
    ],
    warnonly = [:missing_docs, :docs_block],
)

deploydocs(
    repo = "github.com/Statistical-network-analysis-with-Julia/ERGMUserterms.jl.git",
    devbranch = "main",
    versions = [
        "stable" => "dev", # serve dev docs at /stable until a release is tagged
        "dev" => "dev",
    ],
    push_preview = true,
)
