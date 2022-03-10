# archive2dna issues and todo

## Todo

* Acceleration using C complied lib
* More efficient support of line and columns reshaping when appliying Reed Solomon corrections while decoding
* Add directed brute force approach for inner code to compensate for
  frameshift mutations?
* Adapt index genration to mi parameter.
* Optimize memory usage.

## Issues

* Deletion (or to many mutations) in first line leads (in some cases) to decoding error.
* At decoding, make necso and block detection more robust by using the ovearall countdown data

