PBAR - One-line refreshing progress
===================================

:author: M. Fouesneau

Inspired by my usual python version, this progress bar looks like this


.. code::

        description [#---------] k/n  10% [time: 01:23:45s, eta: 01d 00:12:34, 2.7 iters/sec]
 

Example code
------------

 .. code:: c

     size_t n = 10000;
     PBar pb(n);
     pb.set_description("short loop");

     std::cout << "test short loop" << std::endl;
     // neat to also increment the pb during the for declaration!!
     for(size_t i=0; i <= n; ++i, ++pb) {
     	usleep(200);
     }


