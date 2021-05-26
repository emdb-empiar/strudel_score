# The "make" targets are:
#    wheel: build a Python wheel in "dist" directory
#    install: build wheel and install
#    test: run ChimeraX with test commands
#    debug: run ChimeraX with debugging flag set
#    clean: remove files used in building wheel

# Theoretically, you should only need to change
# 1. the CHIMERAX_APP definition (or define in environment),

# Define where ChimeraX is installed.
OS = $(patsubst CYGWIN_NT%,CYGWIN_NT,$(shell uname -s))
# CHIMERAX_APP is the ChimeraX install folder
# We use ?= to let CHIMERAX_APP environment variable override
ifeq ($(OS),Darwin)
# Mac
CHIMERAX_APP ?= /Applications/ChimeraX-1.2.4.app
endif
ifeq ($(OS),CYGWIN_NT)
# Windows
CHIMERAX_APP ?= "c:/Program Files/ChimeraX_daily.app"
endif
# Platform-dependent settings.  Should not need fixing.
# For Windows, we assume Cygwin is being used.
ifeq ($(OS),CYGWIN_NT)
CHIMERAX_EXE = $(CHIMERAX_APP)/bin/ChimeraX.exe
endif
ifeq ($(OS),Darwin)
CHIMERAX_EXE = $(CHIMERAX_APP)/Contents/bin/ChimeraX
PY_EXE = $(CHIMERAX_APP)/Contents/bin/python3.8
PIP_EXE = $(CHIMERAX_APP)/Contents/bin/pip
# $(warning VAR is $(PY_EXE))
endif
ifeq ($(OS),Linux)
CHIMERAX_EXE = $(shell which chimerax)
CHIMERAX_BIN=$(dir $(CHIMERAX_EXE))
PY_EXE = $(CHIMERAX_BIN)python3.7
TMP = $(subst bin/,, $(CHIMERAX_BIN))
PIP_EXE = $(TMP)lib/python3.7/site-packages/pipls

endif
RUN = $(CHIMERAX_EXE) --nogui --exit --cmd

PYSRCS = $(wildcard src/*.py)
CSRCS = $(wildcard src/*.cpp)
SRCS = $(PYSRCS) $(CSRCS)

# If you want to install into ChimeraX.app, add "user false"
# to the "devel build" and "devel install" commands.
# By default, we install for just the current user.

wheel $(WHEEL): bundle_info.xml $(SRCS)
	$(CHIMERAX_EXE) -m PyQt5.pyrcc_main -o src/resources/resources_rc.py src/resources/resources.qrc; \
	$(CHIMERAX_EXE) --nogui --cmd "devel build . ; exit"

install app-install: $(WHEEL)
#	$(PY_EXE) $(PIP_EXE) install pandas;
#	$(PY_EXE) $(PIP_EXE) install biopython;
#	$(PY_EXE) $(PIP_EXE) install mrcfile;
	$(PY_EXE) $(PIP_EXE) install --upgrade $(wildcard lib/*.whl);
	$(CHIMERAX_EXE) -m PyQt5.pyrcc_main -o src/resources/resources_rc.py src/resources/resources.qrc; \
	$(CHIMERAX_EXE) --nogui --cmd "devel install . ; exit"

test:
	$(CHIMERAX_EXE)

debug:
	$(CHIMERAX_EXE) --debug

clean:	clean_generated_files
	if [ -x $(CHIMERAX_EXE) ]; then \
		$(RUN) "devel clean . exit true" ; \
	else \
		rm -rf build dist *.egg-info src/__pycache__ ; \
	fi

pylint:
	$(CHIMERAX_EXE) -m flake8 $(filter %.py, $(SRCS))

# Modify targets below if you have files that need
# to be generated before building the bundle and/or
# removed when cleaning the bundle

generated_files:
	# Generate files

clean_generated_files:
	# Remove generated files
