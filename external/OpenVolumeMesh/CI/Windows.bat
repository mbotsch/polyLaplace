mkdir build-release

cd build-release

IF "%ARCHITECTURE%" == "x64" (
  set ARCH_VS= Win64
  set STRING_ARCH=64-Bit
) else (
  set ARCH_VS=
  set STRING_ARCH=32-Bit
)

IF "%BUILD_PLATFORM%" == "VS2013" (
    set LIBPATH=E:\libs\VS2013
    set GTESTVERSION=gtest-1.6.0
    set GENERATOR=Visual Studio 12%ARCH_VS%
    set VS_PATH="C:\Program Files (x86)\Microsoft Visual Studio 12.0\Common7\IDE\devenv.com"
) 

IF "%BUILD_PLATFORM%" == "VS2015" (
    set LIBPATH=E:\libs\VS2015
    set GTESTVERSION=gtest-1.7.0
    set GENERATOR=Visual Studio 14%ARCH_VS%
    set VS_PATH="C:\Program Files (x86)\Microsoft Visual Studio 14.0\Common7\IDE\devenv.com"

) 

IF "%BUILD_PLATFORM%" == "VS2017" (
    set LIBPATH=E:\libs\VS2017
    set GTESTVERSION=gtest-1.7.0
    set GENERATOR=Visual Studio 15%ARCH_VS%
    set VS_PATH="C:\Program Files (x86)\Microsoft Visual Studio\2017\Professional\Common7\IDE\devenv.com"
)

ECHO "==============================================================="
ECHO "==============================================================="
ECHO "Building with :"
ECHO "ARCHITECTURE        : %ARCHITECTURE%"
ECHO "BUILD_PLATFORM      : %BUILD_PLATFORM%"
ECHO "GTESTVERSION        : %GTESTVERSION%"
ECHO "GENERATOR           : %GENERATOR%"
ECHO "CMAKE_CONFIGURATION : %CMAKE_CONFIGURATION%"
ECHO "VS_PATH             : %VS_PATH%"
ECHO "LIBPATH             : %LIBPATH%"
ECHO "==============================================================="
ECHO "==============================================================="

set GTEST_PREFIX=%LIBPATH%\%ARCHITECTURE%\%GTESTVERSION%
set GTEST_INCLUDE_DIR=%GTEST_PREFIX%\include
set GTEST_LIBRARY=%GTEST_PREFIX%\lib\gtest.lib
set GTEST_MAIN_LIBRARY=%GTEST_PREFIX%\lib\gtest_main.lib


"C:\Program Files\CMake\bin\cmake.exe" -DGTEST_LIBRARY="%GTEST_LIBRARY%" -DGTEST_INCLUDE_DIR="%GTEST_INCLUDE_DIR%" -DGTEST_MAIN_LIBRARY="%GTEST_MAIN_LIBRARY%" -G "%GENERATOR%"  -DCMAKE_BUILD_TYPE=Release %CMAKE_CONFIGURATION% ..

%VS_PATH% /Build "Release" OpenVolumeMesh.sln /Project "ALL_BUILD"

IF %errorlevel% NEQ 0 exit /b %errorlevel%

cd ..

cd src\Unittests\TestFiles

..\..\..\build-release\Unittests\Release\unittests.exe

cd ..\..\..\

IF %errorlevel% NEQ 0 exit /b %errorlevel%


mkdir build-debug

cd build-debug

set GTEST_LIBRARY=%GTEST_PREFIX%\lib\gtestd.lib
set GTEST_MAIN_LIBRARY=%GTEST_PREFIX%\lib\gtest_maind.lib

"C:\Program Files\CMake\bin\cmake.exe" -DGTEST_LIBRARY="%GTEST_LIBRARY%" -DGTEST_INCLUDE_DIR="%GTEST_INCLUDE_DIR%" -DGTEST_MAIN_LIBRARY="%GTEST_MAIN_LIBRARY%" -G "%GENERATOR%"  -DCMAKE_BUILD_TYPE=Debug %CMAKE_CONFIGURATION% ..

%VS_PATH% /Build "Debug" OpenVolumeMesh.sln /Project "ALL_BUILD"

IF %errorlevel% NEQ 0 exit /b %errorlevel%

cd ..

cd src\Unittests\TestFiles

..\..\..\build-debug\Unittests\Debug\unittests.exe


IF %errorlevel% NEQ 0 exit /b %errorlevel%
