SET CHIMERA_PATH=c:\Program Files\ChimeraX
SET CHIMERAX_EXE="%CHIMERA_PATH%\bin\ChimeraX-console.exe"
SET PY_EXE="%CHIMERA_PATH%\bin\python.exe"
SET PYP_EXE="%CHIMERA_PATH%\bin\Scripts\pip3.exe"
SET SITE="%CHIMERA_PATH%\bin\lib\site-packages"
SET SITE="C:\Users\andrei\AppData\Local\UCSF\ChimeraX\1.1\site-packages"
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
    %PY_EXE% %PYP_EXE% install --upgrade -t %SITE% lib\strudel-0.1-py3-none-any.whl
	%CHIMERAX_EXE% -m PyQt5.pyrcc_main -o src/resources/resources_rc.py src/resources/resources.qrc
	%CHIMERAX_EXE% --nogui --cmd "devel install .; exit"
	BREAK
)
)