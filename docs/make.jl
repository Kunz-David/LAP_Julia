push!(LOAD_PATH,"../src/")

using Documenter, LAP_julia, Literate

generated_pages = ["basic_interaction", "registration_functions", "test_registration_functions"]
local_version = get(ENV, "GITHUB_ACTIONS", nothing) == "true"

if local_version == false
    for gen_page_name in generated_pages
       EXAMPLE = joinpath(@__DIR__, "src", "man", "examples", gen_page_name * ".jl")
       OUTPUT = joinpath(@__DIR__, "src/generated")

       function preprocess(str)
           str = replace(str, "example_placeholder" => gen_page_name)
           return str
       end

       Literate.markdown(EXAMPLE, OUTPUT, preprocess = preprocess)
       Literate.notebook(EXAMPLE, OUTPUT, preprocess = preprocess)
       # Literate.script(EXAMPLE, OUTPUT, preprocess = preprocess)
    end
end

makedocs(
    sitename = "Local All Pass registration in Julia",
    pages = if local_version
                ["index_local.md",
                    "Library" => Any[
                        "Public" => "lib/public.md",
                        "Internals" => "lib/private.md"
                    ]
                ]
            else
                ["index.md",
                       "Usage" => Any[
                             map(page_name -> joinpath("generated", page_name * ".md"), generated_pages)...
                          ],
                       "Library" => Any[
                          "Public" => "lib/public.md",
                          "Internals" => "lib/private.md"
                          ]
               ]
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
