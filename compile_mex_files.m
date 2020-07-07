% cd @myAD/private;
mex -v -largeArrayDims COMPFLAGS='-O3 -ftree-vectorize --fopt-info-vec-optimized -fopt-info-missed -Wall' valXder.c;
mex -v -largeArrayDims COMPFLAGS='-O3 -ftree-vectorize -fopt-info-vec-optimized -fopt-info-missed -Wall' matdrivXvecval.c;
% mex -v -largeArrayDims COMPFLAGS='-O2 -ftree-vectorize -ftree-vectorize-verbose=7 -fopt-info-missed -Wall' matvalXmatder.c;
% mex -v -largeArrayDims COMPFLAGS='-O2 -ftree-vectorize -ftree-vectorize-verbose=7 -fopt-info-missed -Wall' dertransp.c;
% cd ../../;
