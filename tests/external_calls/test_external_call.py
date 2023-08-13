"""Test blastn subprocess calls and output parsing."""
from gencdna.external_calls import BinaryExecWithYamlArgs


def test_external_shell_call():
    echo = BinaryExecWithYamlArgs('echo')
    echo_output = echo.run('-n', 'echo-test')
    assert echo_output.stdout.decode('UTF-8') == 'echo-test'
