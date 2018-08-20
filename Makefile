SHELL:=bash
username?=$(shell whoami)

repoRoot:=${REPO_ROOT}

#
# NB: if ${PY_DIRS} is empty we don't want to start picking up stuff in root;
#     hence we feed find what we've asked it to ignore!
#
PY_DIRS = src
pyFiles:=$(shell find _build ${PY_DIRS} -type f -regextype sed -regex ".*\.py" -not -path "*_build/*")


.PHONY: firstTimeSetup
firstTimeSetup:
	sudo apt-get install python3-venv pep8

.PHONY: cleanVenv
cleanVenv:
	rm -rf _build/venv
	python3 -m venv _build/venv

.PHONY: clean
clean:
	rm -rf _build

.PHONY: autopep8
autopep8:
	autopep8 --ignore E203 -i --max-line-length=120 ${pyFiles}

.PHONY: pep8
pep8:
	pep8 --ignore E203,E402 --max-line-length=120 ${pyFiles}

.PHONY: install
install:
	pip install -r requirements.txt