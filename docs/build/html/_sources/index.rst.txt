QESuite
-------

Quantum Espresso Suite (|QESuite|) is a lightweight package to build input files for quantum espresso automatically. Traditionally, 
the inputfiles needed to be written line by line which is prone to errors. |QESuite| offer an intuitive way to create these files
directly within a python interface that will check for the proper variable's names. Moreover, creating this file within a Python script
also enable a simpolification of sequential task such as iteration over different variables, convergence tests, etc.

Additionally, |QESuite| enables you to create input files directly from the symmetries of structures, provides a set of recipes
to perform calculations of specific systems  (2D materials, vdW heterostructures, etc). 

Finally, |QESuite| takes seriously the proper credit to all authors involves in improving the different features, so we provide the 
function :py:meth:`~qe_suite.builder.PWInput.cite`, which identify all functionalities used to build an inputfile and 
returns a citation string.

.. toctree::
   usage
   api
   :maxdepth: 2
   :caption: QESuite Structure:

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`