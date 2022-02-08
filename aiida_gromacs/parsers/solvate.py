# -*- coding: utf-8 -*-
"""
Parsers provided by aiida_gromacs.

Register parsers via the "aiida.parsers" entry point in setup.json.
"""
from aiida.engine import ExitCode
from aiida.parsers.parser import Parser
from aiida.plugins import CalculationFactory
from aiida.common import exceptions
from aiida.orm import SinglefileData

SolvateCalculation = CalculationFactory('gromacs.solvate')


class SolvateParser(Parser):
    """
    Parser class for parsing output of calculation.
    """

    def __init__(self, node):
        """
        Initialize Parser instance

        Checks that the ProcessNode being passed was produced by a SolvateCalculation.

        :param node: ProcessNode of calculation
        :param type node: :class:`aiida.orm.ProcessNode`
        """
        super().__init__(node)
        if not issubclass(node.process_class, SolvateCalculation):
            raise exceptions.ParsingError('Can only parse SolvateCalculation')

    def parse(self, **kwargs):
        """
        Parse outputs, store results in database.

        :returns: an exit code, if parsing fails (or nothing if parsing succeeds)
        """
        outputs = ['stdout', 'outputfile', 'topfile']

        # Check that folder content is as expected
        files_retrieved = self.retrieved.list_object_names()
        files_expected = [self.node.get_option('output_filename'),
                          self.node.inputs.parameters['o'],
                          self.node.inputs.topfile.filename]

        # Note: set(A) <= set(B) checks whether A is a subset of B
        if not set(files_expected) <= set(files_retrieved):
            self.logger.error("Found files '{}', expected to find '{}'".format(
                files_retrieved, files_expected))
            return self.exit_codes.ERROR_MISSING_OUTPUT_FILES

        # add outputs
        for index, thing in enumerate(files_expected):
            self.logger.info("Parsing '{}'".format(thing))
            with self.retrieved.open(thing, 'rb') as handle:
                output_node = SinglefileData(file=handle)
            self.out(outputs[index], output_node)

        return ExitCode(0)
