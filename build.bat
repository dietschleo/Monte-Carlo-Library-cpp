@echo off
cd /d %~dp0
g++ main.cpp functions.cpp -o montecarlo.exe -mconsole
pause
