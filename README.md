# Dev Guide

## Style

For the most part, just try to do what the Fortran best practices say, but
there are a few things I need to suggest for practical reasons.

- No global variable, even in a module. The plan is for the oracles functions
  to be accessible from a higher level language such as Python and Julia. This
  means that everything that a subroutine uses have to be passed as an
  argument.
    

## Shenanigans

In the `si_data`, integers are 64 bits. This works with Intel's toolchain but
not so well with GNU. I tried `fdefault-integer-8`, but I had mixed results. We
can stick with Intel or try to find a way to fix this.
