@echo off
cd /d %~dp0
g++ main.cpp functions.cpp Simulation.cpp RandomNumber.cpp EuropeanOption.cpp CompoundOption.cpp AmericanOption.cpp AsianOption.cpp -o build.exe -mconsole
pause
