import pickle
import random
import shutil
from simulations.common.python.utilities import load_globals, directed_reachable as directed_reachable_python
from simulations.common.python.utilities import d_sep as directed_reachable_R
from pathlib import Path
import numpy as np
import pytest
import subprocess

def generate_random_test_case(n):
    # Generate a random test case with n variables
    variables = list(range(0, n))
    random.shuffle(variables)

    x = variables[0]
    y = variables[1]

    C = random.sample(variables[2:], random.randint(0, n - 3))
    J = random.sample(variables[2:], random.randint(0, n - 3))

    D = np.random.randint(2, size=(n, n))
    B = np.random.randint(2, size=(n, n))
    B = (B + B.T) // 2  # Symmetrize B
    np.fill_diagonal(B, 0)

    U = np.zeros((n, n), dtype=int)

    return {'x': x, 'y': y, 'C': C, 'J': J, 'D': D, 'B': B, 'U': U}


@pytest.mark.parametrize("n", [5])
def test_outputs(n, num_test_cases = 10):
    random.seed(1)
    np.random.seed(0)

    for _ in range(num_test_cases):
        test = generate_random_test_case(n)

        data = {'x': test['x'] + 1,
                'y': test['y'] + 1,
                'C': [c + 1 for c in test['C']],
                'J': [j + 1 for j in test['J']],
                'D': test['D'],
                'B': test['B'],
                'U': test['U']}

        path = Path.cwd()
        p = Path(path, 'data.pkl')
        with p.open('wb') as f:
            pickle.dump(data, f)
        
        # Call R implementation and get the output
        command = f'Rscript simulations/common/R/utilities_main.R {p} > output.txt'
        subprocess.call(command, shell=True)

        r_output_path = Path(path, 'output.txt')
        with r_output_path.open('r') as f:
            r_output = f.read()

        # Call Python implementation and get the output
        python_output = directed_reachable_python(
            test['x'], test['y'], test['C'], test['J'], test['D'], test['B'], test['U']
        )
        
        # Compare the outputs
        assert (r_output[-5:-1] == "TRUE") == python_output

        # Clean up
        p.unlink()
        r_output_path.unlink()
