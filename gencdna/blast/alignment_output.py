"""Interface for running command line blastn executable."""

import io
import subprocess  # noqa: S404

import pandas as pd


class AlignmentOutputParser(object):

    def __init__(
        self,
        column_names: list[str],
        alignment_output: subprocess.CompletedProcess,
    ) -> None:
        self.returncode = alignment_output.returncode
        self.standard_output = alignment_output.stdout.decode('UTF-8')
        self.standard_error = alignment_output.stderr.decode('UTF-8')
        self.column_names = column_names

    def output_as_dataframe(self) -> pd.DataFrame:
        return pd.read_csv(
            io.StringIO(self.standard_output),
            delimiter='\t',
            header=None,
            names=self.column_names,
        )
