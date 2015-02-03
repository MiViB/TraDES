@echo off
if not exist %WINDIR%\ncbi.ini copy ncbi.ini %WINDIR%\ncbi.ini > NUL
rem Make program group and icon now
if "%OS%"=="Windows_NT" goto nticon
rem Windows 95/98 setup
if "%WINDIR%"=="" goto end2
if not exist "%WINDIR%\Start Menu\Programs\TraDES" mkdir "%WINDIR%\Start Menu\Programs\TraDES" > NUL
shortcut -f -t vistraj.exe -n "%WINDIR%\Start Menu\Programs\TraDES\VisTraj" -d . > NUL
goto end
:nticon
if "%ALLUSERSPROFILE%"=="" goto singleprofile
if not exist "%ALLUSERSPROFILE%\Start Menu\Programs\TraDES" mkdir "%ALLUSERSPROFILE%\Start Menu\Programs\TraDES" > NUL
shortcut -f -t vistraj.exe -n "%ALLUSERSPROFILE%\Start Menu\Programs\TraDES\VisTraj" -d . > NUL
goto end
:singleprofile
if "%USERPROFILE%"=="" goto end2
if not exist "%USERPROFILE%\Start Menu\Programs\TraDES" mkdir "%USERPROFILE%\Start Menu\Programs\TraDES" > NUL
shortcut -f -t vistraj.exe -n "%USERPROFILE%\Start Menu\Programs\TraDES\VisTraj" -d . > NUL
:end
echo Program icon successfully added to start menu
:end2
echo Done.  Now type "vistraj" or "foldtraj" at the DOS prompt to begin.
