# Sphinx variable
SPHINXOPTS    =
SPHINXBUILD   = python -msphinx
SPHINXPROJ    = PylProtPredictor
SOURCEDIR     = src/docs
BUILDDIR      = tmp
MINICONDA_URL = https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
MINICONDA	  = $HOME/miniconda3
CONDA_ENV     = PylProtPredictor

ifeq ($(shell uname -s),Darwin)
	MINICONDA_URL=https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
endif

CONDA=$(shell which conda)
ifeq ($(CONDA),)
	CONDA=${HOME}/miniconda3/bin/conda
endif

# Commands
default: help

install-conda: ## install Miniconda
	curl $(MINICONDA_URL) -o miniconda.sh
	bash miniconda.sh -b
.PHONY: install-conda

create-env: ## create conda environment
	if ${CONDA} env list | grep '^${CONDA_ENV}'; then \
	    ${CONDA} env update -f environment.yml; \
	else \
	    ${CONDA} env create -f environment.yml; \
	fi
.PHONY: create-env

ACTIVATE_ENV = source activate ${CONDA_ENV}

init: ## install the requirements
	python setup.py install
.PHONY: init

develop: init ## setup develop mode
	python setup.py develop
.PHONY: develop

test: ## run the tests
	flake8 --exclude=.git,build,.eggs --ignore=E501 .
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
