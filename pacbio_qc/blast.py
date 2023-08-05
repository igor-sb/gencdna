"""Interface for running command line blastn executable."""

import io
import subprocess  # noqa: S404

import pandas as pd
import yaml


def is_executable_available(executable_path: str) -> bool:
    cmd = subprocess.run([executable_path], capture_output=True)  # noqa: S603
    return cmd.returncode == 1


class Blast(object):

    def __init__(
        self,
        config_yml: str = 'config/blast.yml',
        executable_path: str = 'blastn',
    ) -> None:
        with open(config_yml) as config_file:
            self.config = yaml.safe_load(config_file)
        if not is_executable_available(executable_path):
            raise FileNotFoundError(executable_path)
        self.executable_path = executable_path

    def create_blast_command(self) -> list[str]:
        cmd = [self.executable_path]
        for arg_flag, arg_value in self.config['alignment'].items():
            cmd.extend([
                '-{arg_flag}'.format(arg_flag=arg_flag),
                str(arg_value),
            ])
        return cmd

    def run(self, query: str, subject: str) -> subprocess.CompletedProcess:
        cmd = self.create_blast_command()
        cmd.extend(['-query', query, '-subject', subject])
        return subprocess.run(cmd, capture_output=True)  # noqa: S603


class BlastOutputParser(object):

    def __init__(
        self,
        blast_output: subprocess.CompletedProcess,
        config_yml='config/blast.yml',
    ) -> None:
        self.returncode = blast_output.returncode
        self.standard_output = blast_output.stdout.decode('UTF-8')
        self.standard_error = blast_output.stderr.decode('UTF-8')
        with open(config_yml) as config:
            self.column_names = yaml.safe_load(config)['outfmt6_column_names']

    def output_as_dataframe(self) -> pd.DataFrame:
        return pd.read_csv(
            io.StringIO(self.standard_output),
            delimiter='\t',
            header=None,
            names=self.column_names,
        )
