JULIA ?= julia
NPM ?= npm
MARKDOWN = README.md "docs/**/*.md" "examples/**/*.md"

.PHONY: format format-julia format-markdown check-format lint \
	lint-markdown quality

format: format-julia format-markdown

format-julia:
	$(JULIA) --project=dev -e 'using Pkg; Pkg.instantiate()'
	$(JULIA) --project=dev dev/format.jl

format-markdown: node_modules/.package-lock.json
	$(NPM) exec -- prettier --write $(MARKDOWN)

check-format: node_modules/.package-lock.json
	$(JULIA) --project=dev -e 'using Pkg; Pkg.instantiate()'
	$(JULIA) --project=dev dev/check_format.jl
	$(NPM) exec -- prettier --check $(MARKDOWN)

lint:
	$(JULIA) --project=dev -e 'using Pkg; Pkg.instantiate()'
	$(JULIA) --project=dev -m LintApp . \
		--no-progress --verbose --max-warnings 0

lint-markdown: node_modules/.package-lock.json
	$(NPM) exec -- markdownlint-cli2

quality: check-format lint-markdown

node_modules/.package-lock.json: package-lock.json
	$(NPM) ci
