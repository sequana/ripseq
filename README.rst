
.. image:: https://badge.fury.io/py/sequana-ripseq.svg
     :target: https://pypi.python.org/pypi/sequana_ripseq

.. image:: http://joss.theoj.org/papers/10.21105/joss.00352/status.svg
    :target: http://joss.theoj.org/papers/10.21105/joss.00352
    :alt: JOSS (journal of open source software) DOI

.. image:: https://github.com/sequana/ripseq/actions/workflows/main.yml/badge.svg
   :target: https://github.com/sequana/ripseq/actions/workflows    

.. image:: http://joss.theoj.org/papers/10.21105/joss.00352/status.svg
   :target: http://joss.theoj.org/papers/10.21105/joss.00352
   :alt: JOSS (journal of open source software) DOI

.. image:: https://img.shields.io/badge/python-3.8%20%7C%203.9%20%7C3.10-blue.svg
    :target: https://pypi.python.org/pypi/sequana/ripseq
    :alt: Python 3.8 | 3.9 | 3.10



This is is the **ripseq** pipeline from the `Sequana <https://sequana.readthedocs.org>`_ project

:Overview: TODO 
:Input: TODO
:Output: TODO
:Status: draft
:Citation: Cokelaer et al, (2017), ‘Sequana’: a Set of Snakemake NGS pipelines, Journal of Open Source Software, 2(16), 352, JOSS DOI doi:10.21105/joss.00352


Installation
~~~~~~~~~~~~

sequana_fastqc is based on Python3, just install the package as follows::

    pip install sequana_ripseq --upgrade

You will need third-party software such as fastqc. Please see below for details.

Usage
~~~~~

::

    sequana_pipelines_ripseq --help
    sequana_pipelines_ripseq --input-directory DATAPATH 

This creates a directory with the pipeline and configuration file. You will then need 
to execute the pipeline::

    cd ripseq
    sh ripseq.sh  # for a local run

This launch a snakemake pipeline. If you are familiar with snakemake, you can 
retrieve the pipeline itself and its configuration files and then execute the pipeline yourself with specific parameters::

    snakemake -s ripseq.rules -c config.yaml --cores 4 --stats stats.txt

Or use `sequanix <https://sequana.readthedocs.io/en/main/sequanix.html>`_ interface.


Usage with apptainer / singularity::
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

With Apptainer, initiate the working directory as follows::

    sequana_ripseq --use-apptainer

Images are downloaded in the working directory but you can store then in a directory globally (e.g.)::

    sequana_ripseq --use-apptainer --apptainer-prefix ~/.sequana/apptainers

and then::

    cd ripseq
    sh ripseq.sh

if you decide to use snakemake manually, do not forget to add apptainer options::

    snakemake -s ripseq.rules -c config.yaml --cores 4 --stats stats.txt --use-apptainer --apptainer-prefix ~/.sequana/apptainers --apptainer-args "-B /home:/home"

By default, the home is already set for you. Additional binding path can be set using environment variables e.g.::

    export APPTAINER_BINDPATH=" -B /pasteur"

Requirements
~~~~~~~~~~~~

This pipelines requires the following executable(s):

- TODO

.. image:: https://raw.githubusercontent.com/sequana/sequana_ripseq/main/sequana_pipelines/ripseq/dag.png


Details
~~~~~~~~~

This pipeline runs **ripseq** in parallel on the input fastq files (paired or not). 
A brief sequana summary report is also produced.


Rules and configuration details
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here is the `latest documented configuration file <https://raw.githubusercontent.com/sequana/sequana_ripseq/main/sequana_pipelines/ripseq/config.yaml>`_
to be used with the pipeline. Each rule used in the pipeline may have a section in the configuration file. 

Changelog
~~~~~~~~~

========= ====================================================================
Version   Description
========= ====================================================================
0.0.1     **First release.**
========= ====================================================================


Contribute & Code of Conduct
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To contribute to this project, please take a look at the 
`Contributing Guidelines <https://github.com/sequana/sequana/blob/main/CONTRIBUTING.rst>`_ first. Please note that this project is released with a 
`Code of Conduct <https://github.com/sequana/sequana/blob/main/CONDUCT.md>`_. By contributing to this project, you agree to abide by its terms.

