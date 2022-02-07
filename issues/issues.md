# archive2dna issues and todo

## Todo

* Add outer code block management
* Optimize memory usage.
* Add directed brute force approach for inner code to compensate for
  frameshift mutations.
* Adapt index genration to mi parameter.

## Issues

* When last or lasts DNA segments are lost, outer code restauration
  generates a padding of 0 for the last segment. It is removed by the
  internal unzipping at end of decoding process (so it is not a
  problem but unelegant)

