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


class UsearchOutputParser(AlignmentOutputParser):

    def output_as_dataframe(self) -> pd.DataFrame:
        df = pd.read_csv(
            io.StringIO(self.standard_output),
            delimiter='\t',
            header=None,
            names=self.column_names,
        )
        df['query_sequence'] = df.apply(
            lambda row: self._extract_substring(row, 'query'),
            axis=1,
        )
        df['subject_sequence'] = df.apply(
            lambda row: self._extract_substring(row, 'subject'),
            axis=1,
        )
        return df

    def _extract_substring(self, row: pd.Series, col: str):
        start_index: int = row['{col}_start'.format(col=col)] - 1
        stop_index: int = row['{col}_stop'.format(col=col)]
        return row['{col}_sequence'.format(col=col)][start_index:stop_index]
