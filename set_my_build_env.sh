# NOTE: this assumes aarch64 build
export NDK_ROOT="/home/qgao/android-ndk-r17b/"
# aarch64
export ARCH_NAME="aarch64-linux-android"
export ARCH_VER="4.9"
export PATH="$NDK_ROOT/toolchains/${ARCH_NAME}-${ARCH_VER}/prebuilt/linux-x86_64/bin/:$PATH"
echo $PATH | tr ":" "\n"
export SYS_ROOT="$NDK_ROOT/platforms/android-27/arch-arm64/"
export CC="${ARCH_NAME}-gcc --sysroot=$SYS_ROOT"
export LD="${ARCH_NAME}-ld"
export AR="${ARCH_NAME}-ar"
export RANLIB="${ARCH_NAME}-ranlib"
export STRIP="${ARCH_NAME}-strip"
export CFLAGS_INC_DIR="-I$NDK_ROOT/sysroot/usr/include/"
export CFLAGS_INC_DIR="-I$NDK_ROOT/sysroot/usr/include/x86_64-linux-android/ $CFLAGS_INC_DIR"
export CFLAGS_INC_DIR="-I$NDK_ROOT/toolchains/${ARCH_NAME}-${ARCH_VER}/prebuilt/linux-x86_64/lib/gcc/aarch64-linux-android/4.9.x/include/ $CFLAGS_INC_DIR"
 # location to install the lib
export DEST_DIR="$NDK_ROOT/toolchains/${ARCH_NAME}-${ARCH_VER}/prebuilt/linux-x86_64/user/"

