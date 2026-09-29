tar -xjf htslib-1.21.tar.bz2
cd htslib-1.21
autoheader
autoconf
# Use libdeflate when its headers and library are available; otherwise fall back to zlib.
./configure --prefix=`pwd` --with-libdeflate=check
make
make install
