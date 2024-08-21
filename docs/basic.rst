Basic Usage
===========

Required input is a text file of a list of paths to fasta files ``genomes_list`` and a path to an output directory ``output_dir``. All the genomes must be circular and complete. If you have all your genomes in one directory, you can navigate to that directory and generate ``genomes_list`` by running

.. code-block:: console

    ls -d -1 $PWD/*.fasta > input.txt

Then usage is

.. code-block:: console

    pling input.txt output_dir align

