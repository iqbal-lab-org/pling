Basic Usage
===========

Required input is a text file of a list of paths to fasta files ``genomes_list`` and a path to an output directory ``output_dir``. By default, pling assumes the input genomes are circular and complete; use ``--topology`` for mixed circular and linear inputs, or ``--regions`` when clustering regions rather than complete genomes. If you have all your genomes in one directory, you can navigate to that directory and generate ``genomes_list`` by running

.. code-block:: console

    ls -d -1 $PWD/*.fasta > input.txt

Then usage is

.. code-block:: console

    pling cluster align input.txt output_dir
