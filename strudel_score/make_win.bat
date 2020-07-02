SET CHIMERA_PATH=c:\Program Files\ChimeraX
SET CHIMERAX_EXE="%CHIMERA_PATH%\bin\ChimeraX-console.exe"
REM SET PY_EXE="%CHIMERA_PATH%\bin\python.exe"
REM SET PYP_EXE="%CHIMERA_PATH%\bin\Scripts\pip3.exe"
REM SET SITE="%CHIMERA_PATH%\bin\lib\site-packages"

for %%A in (%*) DO (
IF "%%A" == "clean" (
	%CHIMERAX_EXE% --nogui --cmd "devel clean .; exit"
	BREAK
)
)

for %%A in (%*) DO (
IF "%%A" == "app-install" (
	attrib -r +s drive:
REM	%PY_EXE% %PYP_EXE% install -t %SITE% pandas
REM	%PY_EXE% %PYP_EXE% install -t %SITE% biopython
REM	%PY_EXE% %PYP_EXE% install -t %SITE% mrcfile
	%CHIMERAX_EXE% -m PyQt5.pyrcc_main -o src/resources/resources_rc.py src/resources/resources.qrc
	%CHIMERAX_EXE% --nogui --cmd "devel install .; exit"
	BREAK
)
)