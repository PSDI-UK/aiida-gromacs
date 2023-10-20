#!/usr/bin/env python
"""CLI utility to run gmx grompp with AiiDA.

Usage: gmx_grompp --help
"""

import click
import os
import re

from aiida import cmdline, engine
from aiida.plugins import CalculationFactory, DataFactory

from aiida_gromacs import helpers
from aiida_gromacs.utils import searchprevious


def itp_finder(mdpfile, topfile):
    """ Extract included files from the topology file.
    
    This method will grab a list of all includes and then try to sort them
    based on if they are in C directive tags and if flags in the MDP affect
    them.
    """

    # Grab the file contents.
    top = topfile.get_content()
    mdp = mdpfile.get_content()

    # find all include statements in top
    found_include = re.findall(r'#include "([\S\s]*?)"', top, re.DOTALL)

    # find ifdef and define statements in top
    found_ifdef = re.findall(r'#ifdef ([\S\s]*?)\n', top, re.DOTALL)
    defines = re.findall(r'#define ([\S\s]*?)\n', top, re.DOTALL)

    # find defines in mdp file.
    found_def = re.search(r'define([\S\s]*?)\n', mdp, re.DOTALL)
    if (found_def is not None): 
        found_def = re.findall(r'-D([\S\s]*?) ', found_def.group(), re.DOTALL)
        defines.extend(found_def)

    # find ifndef and undefine statements in top
    found_ifndef = re.findall(r'#ifndef ([\S\s]*?)\n', top, re.DOTALL)
    undefines = re.findall(r'#undef ([\S\s]*?)\n', top, re.DOTALL)

    # Extract includes between ifdefs and tag them against vars. 
    for item in found_ifdef:

        # Find the directive
        ifdef_tag = re.search(f'#ifdef {item}([\S\s]*?)#endif', top, re.DOTALL)

        # Find the includes
        if (ifdef_tag is not None) and ((item not in defines) or (item in undefines)):
            ifdef_includes = re.findall(r'#include "([\S\s]*?)"', ifdef_tag.group(), re.DOTALL)
            found_include = list(set(found_include) - set(ifdef_includes))

    # Extract includes between ifndefs and tag them against vars.
    for item in found_ifndef:

        # Find the directive
        ifndef_tag = re.search(f'#ifndef {item}([\S\s]*?)#endif', top, re.DOTALL)

        # Find the includes
        if (ifndef_tag is not None) and ((item in defines) and (item not in undefines)):
            ifndef_includes = re.findall(r'#include "([\S\s]*?)"', ifndef_tag.group(), re.DOTALL)
            found_include = list(set(found_include) - set(ifndef_includes))

    if found_include:
        return gmx_blacklist(found_include)
    else:
        return False


def gmx_blacklist(includes):
    """Remove itp files that are part of gromacs itself."""

    blacklist = [
        "amber03.ff",
        "amber94.ff",
        "oplsaa.ff",
        "gromos53a5.ff",
        "amber96.ff",
        "gromos43a2.ff",
        "amber99sb.ff",
        "amber99.ff",
        "gromos54a7.ff",
        "gromos43a1.ff",
        "amberGS.ff",
        "charmm27.ff",
        "amber99sb-ildn.ff",
        "gromos53a6.ff",
        "gromos45a3.ff"
    ]

    # Remove banned directory substrings
    for banned_ipt in blacklist:
        includes = [ x for x in includes if banned_ipt not in x ]

    return includes


def launch(params):
    """Run grompp.

    Uses helpers to add gromacs on localhost to AiiDA on the fly.
    """

    # Prune unused CLI parameters from dict.
    params = {k:v for k,v in params.items() if v != None}

    # dict to hold our calculation data.
    inputs = {
        "metadata": {
            "description": params.pop("description"),
        },
    }

    # If code is not initialised, then setup.
    if "code" in inputs:
        inputs["code"] = params.pop("code")
    else:
        computer = helpers.get_computer()
        inputs["code"] = helpers.get_code(entry_point="gromacs", computer=computer)

    # Prepare input parameters in AiiDA formats.
    SinglefileData = DataFactory("core.singlefile")
    inputs["mdpfile"] = SinglefileData(file=os.path.join(os.getcwd(), params.pop("f")))
    inputs["grofile"] = SinglefileData(file=os.path.join(os.getcwd(), params.pop("c")))
    inputs["topfile"] = SinglefileData(file=os.path.join(os.getcwd(), params.pop("p")))

    # Find itp files.
    itp_files = itp_finder(inputs["mdpfile"], inputs["topfile"])

    # If we have itps then tag them.
    if itp_files is not False:
        # set correct itpfile path for tests
        if "PYTEST_CURRENT_TEST" in os.environ:
            inputs["itpfiles"] = SinglefileData(file=os.path.join(os.getcwd(), 'tests/input_files', itpfile))
        else:
            inputs["itp_files"] = {}

            # Iterate files to assemble a dict of names and paths.
            for i, itpfile in enumerate(itp_files):
                inputs["itp_files"][f"itpfile{i}"] = SinglefileData(file=os.path.join(os.getcwd(), itpfile))

    if "r" in params:
        inputs["r_file"] = SinglefileData(file=os.path.join(os.getcwd(), params.pop("r")))

    if "rb" in params:
        inputs["rb_file"] = SinglefileData(file=os.path.join(os.getcwd(), params.pop("rb")))

    if "n" in params:
        inputs["n_file"] = SinglefileData(file=os.path.join(os.getcwd(), params.pop("n")))
  
    if "t" in params:
        inputs["t_file"] = SinglefileData(file=os.path.join(os.getcwd(), params.pop("t")))

    if "e" in params:
        inputs["e_file"] = SinglefileData(file=os.path.join(os.getcwd(), params.pop("e")))

    if "qmi" in params:
        inputs["qmi_file"] = SinglefileData(file=os.path.join(os.getcwd(), params.pop("qmi")))

    if "ref" in params:
        inputs["ref_file"] = SinglefileData(file=os.path.join(os.getcwd(), params.pop("ref")))

    GromppParameters = DataFactory("gromacs.grompp")
    inputs["parameters"] = GromppParameters(params)

    # check if inputs are outputs from prev processes
    inputs = searchprevious.get_prev_inputs(inputs, ["grofile", "topfile", "mdpfile"])

    # check if a pytest test is running, if so run rather than submit aiida job
    # Note: in order to submit your calculation to the aiida daemon, do:
    # pylint: disable=unused-variable
    if "PYTEST_CURRENT_TEST" in os.environ:
        future = engine.run(CalculationFactory("gromacs.grompp"), **inputs)
    else:
        future = engine.submit(CalculationFactory("gromacs.grompp"), **inputs)


@click.command()
@cmdline.utils.decorators.with_dbenv()
@cmdline.params.options.CODE()
@click.option("--description", default="record grompp data provenance via the aiida_gromacs plugin", type=str, help="Short metadata description")
# Input file options
@click.option("-f", default="grompp.mdp", type=str, help="Input parameter file")
@click.option("-c", required=True, type=str, help="Input structure file")
@click.option("-r", type=str, help="Structure file: gro g96 pdb brk ent esp tpr")
@click.option("-rb", type=str, help="Structure file: gro g96 pdb brk ent esp tpr")
@click.option("-n", type=str, help="Index file")
@click.option("-p", default="topol.top", type=str, help="Topology file")
@click.option("-t", type=str, help="Full precision trajectory: trr cpt tng")
@click.option("-e", type=str, help="Energy file")
@click.option("-qmi", type=str, help="Input file for QM program")
@click.option("-ref", type=str, help="Full precision trajectory: trr cpt tng")
# Output file options
@click.option("-o", default="conf.gro", type=str, help="Output structure file")
@click.option("-po", type=str, help="grompp input file with MD parameters")
@click.option("-pp", type=str, help="Topology file")
@click.option("-imd", type=str, help="Coordinate file in Gromos-87 format")
# Other parameter options
@click.option("-v", type=str, help="Be loud and noisy")
@click.option("-time", type=str, help="Take frame at or first after this time.")
@click.option("-rmvsbds", type=str, help="Remove constant bonded interactions with virtual sites")
@click.option("-maxwarn", type=str, help="Number of allowed warnings during input processing. Not for normal use and may generate unstable systems")
@click.option("-zero", type=str, help="Set parameters for bonded interactions without defaults to zero instead of generating an error")
@click.option("-renum", type=str, help="Renumber atomtypes and minimize number of atomtypes")
def cli(*args, **kwargs):
    # pylint: disable=unused-argument
    """Run example.

    Example usage:

    $ gmx_grompp --code gmx@localhost -f ions.mdp -c 1AKI_solvated.gro -p 1AKI_topology.top -o 1AKI_ions.tpr

    Alternative (automatically tried to create gmx@localhost code, but requires
    gromacs to be installed and available in your environment path):

    $ gmx_grompp -f ions.mdp -c 1AKI_solvated.gro -p 1AKI_topology.top -o 1AKI_ions.tpr

    Help: $ gmx_grompp --help
    """

    launch(kwargs)


if __name__ == "__main__":
    cli()  # pylint: disable=no-value-for-parameter
