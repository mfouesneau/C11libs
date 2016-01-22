A dead simple option parser
===========================


This is a very simplified parser for command line options
It allows to declare options and access them. In addition, it generates the
usual help message.

What this does not do:

* it does not check that an option used was not declare.
* it does not filter the command line to leave only residuals.


.. code:: cpp

   #include "mfopts.hh"

   int main(int argc, char * argv[])
   {
       Options opt(argv[0]);
       opt.add_option("-h,--help", "display help message");
       opt.add_option("-i,--input", "input file");

       opt.parse_options(argc, argv);

       if (opt.has_option("-h")){
           std::cout<< opt.help();
           exit(0);
       }

       if (opt.has_option("--input")){
           std::cout << opt.get_option<std::string>("--input") << std::endl;
           std::cout << "input = " << opt.get_option<double>("--input") << std::endl;
           std::cout << "input = " << opt.get_option<float>("--input") << std::endl;
       }

       return 0;
  }
