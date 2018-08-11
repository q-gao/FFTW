set exe=%1
adb push %exe% /data/hdr+
adb shell chmod 777 /data/hdr+/%exe%
adb shell /data/hdr+/%exe%