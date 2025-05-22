# creates:  nitrogen.txt, ethane.txt, gold.txt
import os
import subprocess
import sys
import tempfile


def run_script_and_get_output(script):
    """Returns the stdout of executing the code in pythonfile
    as a string."""
    script = os.path.join(os.getcwd(), script)
    with tempfile.TemporaryDirectory() as tempdir:
        return subprocess.check_output(
            [sys.executable, script], cwd=tempdir, encoding='utf8'
        )


# Only save the parts relevant to thermochemistry
output = run_script_and_get_output('nitrogen.py')
output = output[output.find('Enthalpy') :]
with open('nitrogen.txt', 'w') as f:
    f.write(output)

output = run_script_and_get_output('ethane.py')
output = output[output.find('Internal') :]
with open('ethane.txt', 'w') as f:
    f.write(output)

output = run_script_and_get_output('gold.py')
output = output[output.find('Internal') :]
with open('gold.txt', 'w') as f:
    f.write(output)
