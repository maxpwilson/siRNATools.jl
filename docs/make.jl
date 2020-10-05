using Documenter
using siRNATools

makedocs(
    sitename = "siRNATools",
    format = Documenter.HTML(),
    modules = [siRNATools],
    pages = [
        "Home" => "index.md",
        "Installation Guide" => "man/install_guide.md",
        "Informatics" => Any[
            "man/Informatics/batch_index.md"
        ],
        "Modeling" => Any[
        ]
    ]
)

deploydocs(
    repo = "github.com/maxpwilson/siRNATools.jl.git"
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
