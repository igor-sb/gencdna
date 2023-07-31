"""Interface for running command line blastn executable."""

import subprocess  # noqa: S404

import yaml


def is_executable_available(executable_path: str) -> bool:
    cmd = subprocess.run([executable_path], capture_output=True)
    return cmd.returncode == 1


class Blast(object):

    def __init__(
        self,
        query: str,
        subject: str,
        config_yml: str = 'config/blast.yml',
    ) -> None:
        with open(config_yml) as config_file:
            self.config = yaml.safe_load(config_file)
        if not is_executable_available(self.config['executable_path']):
            raise FileNotFoundError(self.config['executable_path'])
        self.query = query
        self.subject = subject

    def create_blast_command(self) -> list[str]:
        cmd = [self.config['executable_path']]
        for arg_flag, arg_value in self.config:
            cmd.extend([
                '-{arg_flag}'.format(arg_flag=arg_flag),
                str(arg_value),
            ])
        cmd.extend(['-query', self.query, '-subject', self.subject])
        return cmd

    def run(self) -> subprocess.CompletedProcess:
        cmd = self.create_blast_command()
        return subprocess.run(cmd, capture_output=True)  # noqa: S603
