Parametric Cut Generation
~~~~~~~~~~~~~~~~~~~~~~~~~

Parametric Cut Generation contains the source code for the python package ``parametricCutGen`` based on the dissertation of Acadia Larsen.

``parametricCutGen`` implements a single row optimal cut selection for Mixed Integer Programs over the domain (and restricted domains) of continuous minimal functions with at most :math:`k` breakpoints. This repository is in an alpha state. 

Current Intents
~~~~~~~~~~~~~~~

 - Illustrate concept of explicit optimal cut selection as a proof of concept for MIP solvers.

 - Reproducibliblity of experimental data using a HPC.

 - Demonstrate use of ``passsagemath`` and ``cutgeneratingfunctionology`` in application; in particular illustrate application of cutting edge mathematics to application of MIPs.

 - Documentation and testing is minimal is intended to serve as a supplement to my dissertation. 

Future Directions
~~~~~~~~~~~~~~~~~

Future directions of this repository include further study of the cut selection problems and developing code to implement cut selection problems algorithms at scale.

Installation and use
~~~~~~~~~~~~~~~~~~~~

Pip install requirements ``cutgeneratingfunctionology``, ``minimalFunctionCache``, ``pyscipopt``, ``pplitepy``, and ``scipy``.

Clone the repository and pip install from source. 


License 
~~~~~~~
The code is released under the GNU General Public License, version 2, or any later version as published by the Free Software Foundation.
