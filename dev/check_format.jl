using JuliaFormatter

formatted = format("."; overwrite = false, verbose = true)
exit(formatted ? 0 : 1)
