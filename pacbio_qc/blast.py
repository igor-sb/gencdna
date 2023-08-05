"""Interface for running command line blastn executable."""

import io
import subprocess  # noqa: S404

import pandas as pd
import yaml


class BlastOutputParser(object):

    def __init__(
        self,
        blast_output: subprocess.CompletedProcess,
        config_yml='config/blast_format.yml',
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
