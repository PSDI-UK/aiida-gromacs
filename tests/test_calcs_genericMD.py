""" Test for genericMD calculation

"""


from aiida_gromacs.utils import searchprevious


def test_process(gromacs_code):
    """Test running a genericMD calculation using the pdb2gmx command as
    an example.
    Note: this does not test that the expected outputs are created of
    output parsing"""

    # pylint: disable=unused-variable
    result, output_dir = searchprevious.run_genericMD_pdb2gmx(gromacs_code)

    assert "pdb2gmx_1AKI_forcefield_gro" in result
    assert "pdb2gmx_1AKI_topology_top" in result
    assert "pdb2gmx_1AKI_restraints_itp" in result


def test_file_name_match(gromacs_code):
    """Test that the file names returned match what was specified on inputs."""

    # pylint: disable=unused-variable
    result, output_dir = searchprevious.run_genericMD_pdb2gmx(gromacs_code)

    assert (
        result["pdb2gmx_1AKI_forcefield_gro"].list_object_names()[0]
        == "pdb2gmx_1AKI_forcefield.gro"
    )
    assert (
        result["pdb2gmx_1AKI_topology_top"].list_object_names()[0]
        == "pdb2gmx_1AKI_topology.top"
    )
    assert (
        result["pdb2gmx_1AKI_restraints_itp"].list_object_names()[0]
        == "pdb2gmx_1AKI_restraints.itp"
    )
