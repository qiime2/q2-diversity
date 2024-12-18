.PHONY: all lint test test-cov viz-alpha-correlation viz-alpha-group-significance viz-alpha-rarefaction install dev clean distclean

PYTHON ?= python

all: viz-alpha-correlation viz-alpha-group-significance viz-alpha-rarefaction

lint:
	q2lint
	flake8

test: all
	py.test

test-cov: all
	py.test --cov=q2_diversity

q2_diversity/_alpha/alpha_correlation_assets/dist:
	cd q2_diversity/_alpha/alpha_correlation_assets && \
	npm install --no-save && \
	npm run build && \
	cp licenses/* dist/

q2_diversity/_alpha/alpha_group_significance_assets/dist:
	cd q2_diversity/_alpha/alpha_group_significance_assets && \
	npm install --no-save && \
	npm run build && \
	cp licenses/* dist/

q2_diversity/_alpha/alpha_rarefaction_assets/dist:
	cd q2_diversity/_alpha/alpha_rarefaction_assets && \
	npm install --no-save && \
	npm run build && \
	cp licenses/* dist/

viz-alpha-correlation: q2_diversity/_alpha/alpha_correlation_assets/dist
viz-alpha-group-significance: q2_diversity/_alpha/alpha_group_significance_assets/dist
viz-alpha-rarefaction: q2_diversity/_alpha/alpha_rarefaction_assets/dist

install: all
	$(PYTHON) -m pip install -v .

dev: all
	pip install -e .

clean: distclean
	rm -rf q2_diversity/_alpha/alpha_correlation_assets/node_modules
	rm -rf q2_diversity/_alpha/alpha_group_significance_assets/node_modules
	rm -rf q2_diversity/_alpha/alpha_rarefaction_assets/node_modules

distclean:
	rm -rf q2_diversity/_alpha/alpha_correlation_assets/dist
	rm -rf q2_diversity/_alpha/alpha_group_significance_assets/dist
	rm -rf q2_diversity/_alpha/alpha_rarefaction_assets/dist
