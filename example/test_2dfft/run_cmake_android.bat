set CMAKE=C:\Apps\Android\sdk\cmake\3.6.4111459\bin\cmake.exe
set  MAKE=C:\Apps\Android\sdk\cmake\3.6.4111459\bin\ninja.exe
set ANDROID_NDK_HOME=C:\Apps\Android\sdk\ndk-bundle
REM need to install cmake from Android SDK Manager
set CMAKE_TC=%ANDROID_NDK_HOME%\build\cmake\android.toolchain.cmake

%CMAKE% -DCMAKE_BUILD_TYPE=Release^
        -DCMAKE_TOOLCHAIN_FILE=%CMAKE_TC%^
        -DANDROID_ABI=arm64-v8a^
        -DCMAKE_MAKE_PROGRAM=%MAKE%^
        -G Ninja^
        -DDEBUG_LEVEL_INFO_LOW=1^
        ..

echo Use ninja to make the build
echo "C:\Apps\Android\sdk\cmake\3.6.4111459\bin\ninja.exe"