========================
Barron Optical Flow Code
========================

`Original code <http://www.csd.uwo.ca/faculty/barron/FTP/>`_ modified to run with modern compilers.

`Examples <https://scivision.co/barron1994opticalflow/>`_

Build
=====

Prereq
------

* Linux: ``apt install gcc gfortran libopenblas-dev``
* Mac: ``brew install gcc openblas``

::

    cd bin
    cmake ..
    make


Horn-Schunk
===========
::

    ./horn new2binarytreet. 0.5 1.5 20 100 ../TESTDATA/TREE_DATA/TRANS out/ -B 150 150 -MH -T 5.0 
    
    ./flow2ps out/horn.modified.new2binarytreet.F-5.00 out/hornOF.ps

    evince out/hornOF.ps
    
Lucas-Kanade
============
::

  ./lucas new2binarytreet. 1.5 20 1.0 ../TESTDATA/TREE_DATA/TRANS out/ -B 150 150
  
  ./flow2ps out/lucas.new2binarytreet.20F-1.00-1.5 out/lucas20F.ps
  
  evince out/lucas20F.ps
