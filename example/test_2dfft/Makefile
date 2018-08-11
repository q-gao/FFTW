PRJ=simple_example
CC=aarch64-linux-android-gcc
NDK_ROOT="/home/qgao/android-ndk-r17b"
# aarch64
ARCH_NAME="aarch64-linux-android"
ARCH_VER="4.9"
SYS_ROOT="${NDK_ROOT}/platforms/android-27/arch-arm64/"
CC="${ARCH_NAME}-gcc --sysroot=${SYS_ROOT}"
LD="${ARCH_NAME}-ld"
AR="${ARCH_NAME}-ar"
RANLIB="${ARCH_NAME}-ranlib"
STRIP="${ARCH_NAME}-strip"
USER_LIB_ROOT=/home/qgao/android-ndk-r17b//toolchains/aarch64-linux-android-4.9/prebuilt/linux-x86_64/user/
CFLAGS_INC_DIR="-I${NDK_ROOT}/sysroot/usr/include/"
CFLAGS_INC_DIR+=" -I${NDK_ROOT}/sysroot/usr/include/x86_64-linux-android/"
CFLAGS_INC_DIR+=" -I${NDK_ROOT}/toolchains/${ARCH_NAME}-${ARCH_VER}/prebuilt/linux-x86_64/lib/gcc/aarch64-linux-android/4.9.x/include/"
CFLAGS_INC_DIR+=" -I${USER_LIB_ROOT}/include"

${PRJ}: ${PRJ}.c
	${CC} -O3 $CFLAGS_INC_DIR $<  -L${USER_LIB_ROOT}/lib -lfftw3 -o $@ 

run: ./${PRJ}
	./${PRJ}

clean:
	rm -f *~ ${PRJ}.o

clean_all: clean
	rm -f ${PRJ}

emacs:
	emacs -nw Makefile ${PRJ}.c 

tar: clean_all
	tar -cvf ../${PRJ}.tar ../${PRJ}

