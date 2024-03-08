# Dev Guide


## Shenanigans

MKL's adds a `-fdefault-integer-8` to the build option, which crashes when
reading from binary files with `integer 4` and vice versa. I decided to use
`-fdefault-integer-8` because I don't know how to unset it when using MKL.

