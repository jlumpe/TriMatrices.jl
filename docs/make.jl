using Documenter, TriMatrices


makedocs(;
	modules=[TriMatrices],
	authors="Jared Lumpe",
	repo="https://github.com/jlumpe/TriMatrices.jl/blob/{commit}{path}#L{line}",
	sitename="TriMatrices.jl",
	pages=[
		"index.md",
		"api.md",
	],
)
