#!/bin/sh
set -eu

cd "$(dirname "$0")"
tar -xjf htslib-1.24.tar.bz2
cd htslib-1.24

# Match CMake's default macOS architecture, including on Apple Silicon.
if [ "$(uname -s)" = Darwin ]; then
    HTSLIB_ARCH=${HTSLIB_ARCH:-x86_64}
    # Some configure probes replace CFLAGS, so the architecture belongs in CC.
    CC="${CC:-cc} -arch $HTSLIB_ARCH"
    CFLAGS="${CFLAGS:--O2}"
    LDFLAGS="${LDFLAGS:-} -arch $HTSLIB_ARCH"
    if [ "$HTSLIB_ARCH" = x86_64 ]; then
        CPPFLAGS="${CPPFLAGS:-} -I/usr/local/include"
        LDFLAGS="$LDFLAGS -L/usr/local/lib"
        export CPPFLAGS
    fi
    export CC CFLAGS LDFLAGS
fi

# Use libdeflate when its headers and library are available; otherwise fall back to zlib.
./configure --prefix="$PWD" --with-libdeflate=check "$@"
make -j "${JOBS:-2}"
make install
