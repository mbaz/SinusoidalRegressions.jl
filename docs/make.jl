# Workaround for GR warnings
ENV["GKSwstype"] = "100"

using Documenter, SinusoidalRegressions

# suitable for local browsing without a web server
makedocs(sitename="SinusoidalRegressions.jl",
         build = "stable", # change for other versions
         format = Documenter.HTML(prettyurls = false)
        )

# suitable for hosting in github
makedocs(sitename="SinusoidalRegressions.jl",
         build = "stable" # change for other versions
        )
