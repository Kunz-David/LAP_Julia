push!(LOAD_PATH,"../src/")

using Documenter, LAP_julia

makedocs(
    sitename = "Local All Pass registration in Julia",
    pages = ["index.md",
             "Usage" => Any[
                "Examples" => Any[
                   "man/examples/examples.md"
                   ]
                ],
             "Library" => Any[
                "Public" => "lib/public.md",
                "Internals" => "lib/private.md"
                ]
             ],
    # see here https://juliadocs.github.io/Documenter.jl/stable/man/guide/
    format = Documenter.HTML(prettyurls = get(ENV, "GITHUB_ACTIONS", nothing) == "true"),
    authors = "David Kunz"
    )


deploydocs(
    repo = "github.com/Kunz-David/LAP_Julia.jl.git",
)















#     modules = [LAP_julia],
#     clean = true,
#     doctest = false,
#     checkdocs = :all,
#     format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
#     authors = "Gregory L. Wagner and Navid C. Constantinou",
#     sitename = "LAP_julia.jl",
#
#     pages = Any["Home" => "index.md",
#                 "Code Basics" => "basics.md",
#                 "Forcing" => "forcing.md",
#                 "DocStrings" => Any["man/types.md", "man/functions.md"]
#                  ]
# )
#
# deploydocs(
#     repo = "github.com/LAP_julia/LAP_julia.jl.git",
# )
