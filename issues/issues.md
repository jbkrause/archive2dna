# archive2dna issues and todo

## Todo

* Add outer code block management
* Remove dependency to numpy.
* Optimize memory usage.
* Add directed brute force approach for inner code to compensate for
  frameshift mutations.

## Issues

* When last or lasts DNA segments are lost, outer code restauration
  generates a padding of 0 for the last segment. Not a problem for zip
  files but inelegant.
