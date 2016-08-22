MFIO - Trying to make ascii IO simpler.
=======================================

`AsciiReader` and `AsciiWriter` make very simple proxies to reading and writing
Ascii data. They both also handle `stdin` and `stdout`.

@author M. Fouesneau

example
-------

.. code:: c

         // open data file 
         AsciiReader rd("file.csv");
         AsciiWriter wt("-");        // will write on stdout

         // set delimiter and comment char;
         rd.set_delimiter(",");
         rd.set_comment("#");

         // wt.set_delimiter(" "); // space is default value

         // select columns to use
         rd.set_columns({0, 1, 3});

         // define data field values and types
         double ra;
         double dec;
         double gmag;

         // while rd can read lines from the input file
         while(++rd){
                rd.get_next_value<double>(ra);
                rd.get_next_value<double>(dec);
                rd.get_next_value<double>(gmag);
                wt << {ra, dec, gmag};  // write a vector in one command
                ++ wt;                  // new line
         }

Complex string formatting exports can be done through `io::string_format`, which
gives the same API as `printf`.


AsciiReader
-----------

Reads in ascii data.

Notes
~~~~~

* `count_lines` is based on GNU wc implementation for speed. It can skip
  commented lines and be used to prepare memory arrays. It also caches the
  results after the first call.

* reading a new line (`operator ++`) will always skip commented lines.
 
* parsing values through `get_next_value` will ignore non-used columns (if set).


+--------------------------+---------------------------------------------------+
| `AsciiReader`            | Constructor                                       |
+--------------------------+---------------------------------------------------+
| `count_lines`            | count the number of lines ignoring commented ones |
+--------------------------+---------------------------------------------------+
| `get_columns`            | returns which columns are used                    |
+--------------------------+---------------------------------------------------+
| `get_comment`            | which comment character                           |
+--------------------------+---------------------------------------------------+
| `get_delimiter`          | returns the delimiter                             |
+--------------------------+---------------------------------------------------+
| `get_header`             | returns the first non-commented line              |
+--------------------------+---------------------------------------------------+
| `get_next_value`         | parse next value into variable                    |
+--------------------------+---------------------------------------------------+
| `get_number_of_fields`   | counts the number of fields                       |
+--------------------------+---------------------------------------------------+
| `operator ++`            | buffer the next line                              |
+--------------------------+---------------------------------------------------+
| `set_columns`            | set columns to use                                |
+--------------------------+---------------------------------------------------+
| `set_comment`            | set comment definition                            |
+--------------------------+---------------------------------------------------+
| `set_delimiter`          | set the delimiter                                 | 
+--------------------------+---------------------------------------------------+


AsciiWriter
----------

+-----------------+-------------------------------+
| `AsciiWriter`   | constructor                   | 
+-----------------+-------------------------------+
| `get_comment`   | which comment character       | 
+-----------------+-------------------------------+
| `get_delimiter` | returns the delimiter         |
+-----------------+-------------------------------+
| `operator ++`   | start a new line              |
+-----------------+-------------------------------+
| `operator <<`   | input data (value or vectors) |
+-----------------+-------------------------------+
| `set_comment`   | set comment definition        |
+-----------------+-------------------------------+
| `set_delimiter` | set the delimiter             | 
+-----------------+-------------------------------+
