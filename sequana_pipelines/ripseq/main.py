#
#  This file is part of Sequana software
#
#  Copyright (c) 2016-2021 - Sequana Dev Team (https://sequana.readthedocs.io)
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  Website:       https://github.com/sequana/sequana
#  Documentation: http://sequana.readthedocs.io
#  Contributors:  https://github.com/sequana/sequana/graphs/contributors
##############################################################################
import os

import rich_click as click
import click_completion

click_completion.init()


from sequana_pipetools.options import *
from sequana_pipetools import SequanaManager


NAME = "ripseq"

help = init_click(NAME, groups={
    "Pipeline Specific": [
        "--rmethod", "--skip-multiqc"],
        }
)



@click.command(context_settings=help)
@include_options_from(ClickSnakemakeOptions, working_directory=NAME)
@include_options_from(ClickSlurmOptions)
@include_options_from(ClickInputOptions, add_input_readtag=True)
@include_options_from(ClickGeneralOptions)
@click.option(
    "--reference",
    "reference",
    default="TODO",
    type=click.Path(),
    show_default=True,
    help="""TODO""",
)
@click.option(
    "--todo",
    "todo",
    default="TODO",
    type=click.Choice(["TODO", "OTHERS"]),
    show_default=True,
    help="""TODO""",
)
@click.option(
    "--skip-todo",
    is_flag=False,
    help="""a flag""",
)
def main(**options):


    # the real stuff is here
    manager = SequanaManager(options, NAME)
    options = manager.options

    # creates the working directory
    manager.setup()

    # Fill the config file with data and specific options
    cfg = manager.config.config
    cfg.input_pattern = options.input_pattern
    cfg.input_directory = os.path.abspath(options.input_directory)
    #cfg.general.todo = options.todo

    cfg.reference = options.reference

    manager.exists(cfg.input_directory)

    # finalise the command and save it; copy the snakemake. update the config
    # file and save it.
    manager.teardown()


if __name__ == "__main__":
    main()
