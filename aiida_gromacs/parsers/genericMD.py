"""
Parsers provided by aiida_gromacs.

This parser saves outputted files from a generic command.
"""

import os
import re

from aiida.common import exceptions
from aiida.engine import ExitCode
from aiida.orm import SinglefileData
from aiida.parsers.parser import Parser
from aiida.plugins import CalculationFactory

# entry point string under which the parser class is registered:
GenericCalculation = CalculationFactory("genericMD")


class GenericParser(Parser):
    """
    Parser class for parsing output of genericMD calculation from which
    the retrieved outputs files from the calcjob and the nodes of finished
    calculation can be accessed.
    """

    def __init__(self, node):
        """
        Initialize Parser instance

        Checks that the ProcessNode being passed was produced by a 
        GenericCalculation.

        :param node: ProcessNode of calculation
        :param type node: :class:`aiida.orm.nodes.process.process.ProcessNode`
        """
        super().__init__(node)
        if not issubclass(node.process_class, GenericCalculation):
            raise exceptions.ParsingError("Can only parse GenericCalculation")

    def parse(self, **kwargs):
        """
        Parse outputs, store results in the AiiDA database.

        :returns: an exit code, if parsing fails or the user defined 
            output files are not returned
        """

        # get_option() convenience method is used to get the filename of
        # the output file
        output_filename = self.node.get_option("output_filename")
        # the directory for storing parsed output files
        output_dir = self.node.get_option("output_dir")

        # Check that folder content is as expected
        files_retrieved = self.retrieved.list_object_names()
        files_expected = []  # [output_filename]
        if "output_files" in self.node.inputs:
            for name in self.node.inputs.output_files:
                files_expected.extend([str(name)])

        # Check all outputted files produced have been previously 
        # defined by the user
        for file in files_expected:
            if file not in files_retrieved:
                self.logger.error(
                    f"User defined output file '{file}' not in "
                    f"list of retrieved files '{files_retrieved}'"
                )
                return self.exit_codes.ERROR_UNTRACKED_OUTPUT_FILES

        # passing along the std output file as a SinglefileData node.
        self.logger.info(f"Parsing '{output_filename}'")
        with self.retrieved.open(output_filename, "rb") as handle:
            output_node = SinglefileData(file=handle)

        # return stdout file 
        self.out("log", output_node)

        # passing along all expected output file as SinglefileData nodes.
        for thing in files_expected:
            self.logger.info(f"Parsing '{thing}'")
            with self.retrieved.open(thing, "rb") as handle:
                output_node = SinglefileData(file=handle, filename=thing)
            self.out(self.format_link_label(thing), output_node)

        # parse retrieved files and write them to where command was run
        for thing in files_retrieved:
            self.logger.info(f"Parsing '{thing}'")
            file_path = os.path.join(output_dir, thing)
            # file_path3 = os.path.join(output_dir, f'{thing}-test2.txt')
            try:
                with self.retrieved.open(thing, "rb") as handle:
                    with open(file_path, "wb") as f_out:
                        while True:
                            chunk = handle.read(1024)
                            if not chunk:
                                break
                            f_out.write(chunk)

                # not used yet.
                # test below for parsing log file and saving output as dict
                '''if re.search('.log$', thing):
                    start_string = 'Input Parameters:' #'A ?V ?E ?R ?A ?G ?E ?S'
                    end_string = 'compressibility' #'M ?E ?G ?A ?- ?F ?L ?O ?P ?S'
                    with self.retrieved.open(thing, "r") as file:
                        file_content = file.read()
                        pattern = rf"{start_string}(.*?){end_string}"
                        matches = re.findall(pattern, file_content, re.DOTALL)
                    file_path3 = os.path.join(output_dir, f'{thing}-matched.txt')
                    # matched_text = self._parse_gromacs_top(thing)
                    with open(file_path3, "w") as f_out:
                        f_out.write('matched text:')
                        for match in matches:
                            lines = match.splitlines()
                            for line in lines:
                                f_out.write(f'{line}\n')
                    with open(file_path3, "rb") as f_out:
                        # output_node = SinglefileData(file=f_out, filename=f'{thing}-matched.txt')
                        # output_node = Dict(file={"test": "test"}, filename='dict.txt')
                        output_node = Dict({"test": "test"})
                        output_node.label = 'test-dict'
                    self.out(self.format_link_label('test-dict'), output_node)'''

            except UnicodeDecodeError:
                with self.retrieved.open(thing, "r") as handle:
                    with open(file_path, "w", encoding="utf-8") as f_out:
                        for line in handle.read():
                            f_out.write(line)

        return ExitCode(0)
    
    def _parse_gromacs_top(self, file_path):
        """Not used yet, test for parsing the gromacs tpr file.
        :param file_path: The path and name of gtomacs .log file
        :returns: The required text from the parsed file
        """

        def _find_text_between_strings(file_path, start_string, end_string):
            with self.retrieved.open(file_path, "r") as file:
                file_content = file.read()
            # use re.escape to escape special characters in the strings
            # pattern = rf"{re.escape(start_string)}(.*?){re.escape(end_string)}"
            pattern = rf"{start_string}(.*?){end_string}"
            # use re.DOTALL to make the dot character match newline characters as well
            matches = re.findall(pattern, file_content, re.DOTALL)
            return matches

        #file_path = '1AKI_production.log'
        start_string = 'A ?V ?E ?R ?A ?G ?E ?S'
        end_string = 'M ?E ?G ?A ?- ?F ?L ?O ?P ?S'
        # start_string = '<======  ###############  ==>\n\t<====  A V E R A G E S  ====>\n\t<==  ###############  ======>\n\n' #'A V E R A G E S'
        # end_string = 'M E G A - F L O P S   A C C O U N T I N G'
        matched_text = _find_text_between_strings(file_path, start_string, end_string)
        return matched_text
    


    @staticmethod
    def format_link_label(filename: str) -> str:
        """
        Modified from: https://github.com/sphuber/aiida-shell/blob/master/src/aiida_shell/parsers/shell.py
        Format the link label from a given filename with prefix.
        Valid link labels can only contain alphanumeric characters and
        underscores, without consecutive underscores. So all characters
        that are not alphanumeric or an underscore are converted to
        underscores, where consecutive underscores are merged into one.
        Additional: Label cannot start with a number or underscore.

        :param filename: The filename.
        :returns: The link label.
        """
        filename = filename.split("/")[-1]
        if filename[0].isdigit() or filename[0] == "_":
            filename = "output_files_" + filename
        alphanumeric = re.sub("[^0-9a-zA-Z_]+", "_", filename)
        link_label = re.sub("_[_]+", "_", alphanumeric)

        return link_label

