using Documenter, cereal

makedocs(
            sitename    =   "cereal.jl"
            authors     = "Justin C. Feng and Filip Hejda and Sante Carloni",
            pages = [
                "Home"              => "index.md",
                "Evaluation"        => "evaluation.md",
                "Locator Functions" => "locatorfuncs.md"
            ]
        )
