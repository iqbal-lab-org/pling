# Inspired by https://github.com/snakemake/snakefmt/blob/master/Makefile
PROJECT = pling
OS := $(shell uname -s)
VERSION := $(shell poetry version -s)
BOLD := $(shell tput bold)
NORMAL := $(shell tput sgr0)

.PHONY: all
all: install

.PHONY: install
install:
	poetry install

.PHONY: install-ci
install-ci:
	poetry install --no-interaction
	poetry run pling --help

.PHONY: test
test:
	poetry run python -m unittest discover -s tests -t .

.PHONY: build
build:
	poetry build
