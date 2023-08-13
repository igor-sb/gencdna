"""Wrapper for binary executables."""
import subprocess  # noqa: S404

import yaml


class BinaryExecutable(object):

    def __init__(self, executable: str, config_yml: str = '') -> None:
        if config_yml == '':
            self.config: dict[str, dict] = {'arguments': {}}
        else:
            with open(config_yml) as config_file:
                self.config = yaml.safe_load(config_file)
        self.executable = executable

    def create_command(self) -> list[str]:
        cmd = [self.executable]
        for arg_flag, arg_value in self.config['arguments'].items():
            cmd.append('{arg_flag}'.format(arg_flag=arg_flag))
            if str(arg_value):
                cmd.append(str(arg_value))
        return cmd

    def run(self, *args) -> subprocess.CompletedProcess:
        cmd = self.create_command()
        cmd.extend(args)
        return subprocess.run(cmd, capture_output=True)  # noqa: S603
