using Documenter, cereal

makedocs(
    sitename    =   "cereal.jl",
    format      = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    authors     = "Justin C. Feng and Filip Hejda and Sante Carloni",
    pages = [
                "Home"              => "index.md",
                "Locator Functions" => "locatorfuncs.md",
                "Evaluation"        => "evaluation.md",
            ]
        )
