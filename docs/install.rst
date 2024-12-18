Installation
============

The sequenza-utils code is hosted in `BitBucket`_.

Supported Python version are 2.7+, including python3, and pypy.

Sequenza-utils can be installed either via the Python Package Index (PyPI)
or from the git repository.

The package uses the external command line tools `samtools`_ and `tabix`_.
For the package to function correctly such programs need to be installed in the system

Latest release via PyPI
-----------------------

To install the latest release via PyPI using pip::

    pip install sequenza-utils


Development version
-------------------
To use the test suite for the package is necessary to install also `bwa`_
Using the latest development version directly from the git repository::

    git clone https://bitbucket.org/sequenzatools/sequenza-utils
    cd sequenza-utils
    python setup.py test
    python setup.py install

.. _`BitBucket`: https://bitbucket.org/sequenzatools/sequenza-utils
.. _`samtools`: http://samtools.sourceforge.net
.. _`tabix`: http://www.htslib.org/doc/tabix.html
.. _`bwa`: http://bio-bwa.sourceforge.net
