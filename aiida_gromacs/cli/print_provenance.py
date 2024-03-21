#!/usr/bin/env python
"""
Show text output of the provenance of commands performed in gmx
"""
import click
from aiida_gromacs.utils.displayprovenance import show_provenance_text

def launch(params):
    """Run print_provenance
    """
    show_provenance_text()

@click.command()
def cli(*args, **kwargs):
    # pylint: disable=unused-argument
    """Print out the provenance of aiida-gromacs processes run on current
    loaded aiida profile

    Example usage:

    $ print_provenance

    Help: $ print_provenance --help
    """
    launch(kwargs)
    

if __name__ == "__main__":
    cli()  # pylint: disable=no-value-for-parameter