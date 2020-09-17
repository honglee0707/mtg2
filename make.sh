
  #yum install glibc-static
  #yum install zlib-devel
  #yum install zlib-static

  source /opt/intel/bin/compilervars.sh intel64
  source /opt/intel/mkl/bin/mklvars.sh intel64

  # set unlimited stack size
  ulimit -s unlimited

  rm *.o
  make -f Makefile-mtg2

