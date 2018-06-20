# Sphinx variable
SPHINXOPTS    =
SPHINXBUILD   = python -msphinx
SPHINXPROJ    = PylProtPredictor
SOURCEDIR     = src/docs
BUILDDIR      = tmp
MINICONDA_URL = https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
MINICONDA	  = $HOME/miniconda3

ifeq ($(shell uname -s),Darwin)
	MINICONDA_URL=https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
endif

# Commands
default: help

init: ## install the requirements
	python setup.py install
.PHONY: init

develop: init ## setup develop mode
	python setup.py develop
.PHONY: develop

install-conda: ## install Miniconda
	wget $(MINICONDA_URL) -O miniconda.sh
	bash miniconda.sh -b
	conda config --set show_channel_urls yes --set always_yes yes
	conda update conda conda-env
.PHONY: install-conda

configure-conda: ## configure conda channels
	conda config --system --add channels conda-forge
	conda config --system --add channels defaults
	conda config --system --add channels r
	conda config --system --add channels bioconda
.PHONY: configure-conda

create-env: ## create conda environment
	conda env create --name PylProtPredictor --file environment.yml
.PHONY: create-env

test: ## run the tests
	flake8 --exclude=.git,build --ignore=E501 .
	pytest --cov=pylprotpredictor tests/
.PHONY: test

doc: ## generate HTML documentation
	@$(SPHINXBUILD) -M html "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
	rm -rf docs
	mv "$(BUILDDIR)/html" docs
	rm -rf docs/_sources
	rm -rf tmp
.PHONY: doc

coverage: ## send test coverage to Code coverage
	codecov -t e0873aae-ab97-4e87-a000-1bf66c75b4ee
.PHONY: coverage

help:
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'
.PHONY: help
