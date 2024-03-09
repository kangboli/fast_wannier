# Dev Guide

## Style

For the most part, just try to do what the Fortran best practices say, but
there are a few things I need to suggest for practical reasons.

- No global variable, even in a module. The plan is for the oracles functions
  to be accessible from a higher level language such as Python and Julia. This
  means that everything that a subroutine uses have to be passed as an
  argument.
    

## Shenanigans

The data I included are Fortran binary files. The integers are saved as 64 bits
and they have to match the default integer length of the compiler. `gfortran`'s
default is 32 bits whereas `ifx`'s default is 32 bits.

Fortran silently deallocates array whose intent is out. This prevents arrays
from being allocated and deallocated in the same scope and paves the road to
Helheim. For now, I'm using `load` for subroutines that can silently deallocate
memory, but I don't know what is the right way to get around this.
