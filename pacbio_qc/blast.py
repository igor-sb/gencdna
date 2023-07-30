"""Interface for running command line blastn executable."""

import os
import subprocess  # noqa: S404

import yaml


class Blast(object):

    def __init__(
        self,
        executable_path: str,
        config_yml: str = 'config/blast.yml',
    ) -> None:
        with open(config_yml) as config_file:
            self.config = yaml.safe_load(config_file)
        if not os.path.exists(executable_path):
            raise FileNotFoundError(executable_path)
        self.executable_path = executable_path

    def create_command(self) -> list[str]:
        cmd = [self.executable_path]
        for arg_flag, arg_value in self.config:
            cmd.extend([
                '-{arg_flag}'.format(arg_flag=arg_flag),
                str(arg_value),
            ])
        return cmd

    def run(self) -> str: # check if its a str or if even we want to return anything
        cmd = self.create_command()
        blast_result = subprocess.run(cmd, capture_output=True)  # noqa: S603
        return blast_result.stdout
