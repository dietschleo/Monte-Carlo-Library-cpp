@echo off
cd /d %~dp0
g++ main.cpp functions.cpp RandomNumber.cpp -o build.exe -mconsole
pause
