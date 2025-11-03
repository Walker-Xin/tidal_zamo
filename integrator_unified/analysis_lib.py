import numpy as np
import pandas as pd
from typing import Dict, Any, Optional


def check_stable_r(data, threshold=0.1, hysterisis=0.1):
    """
    Check for a period of stable r using a moving frame.
    Frame size is set to 500 by default.
    """

    diffs = np.diff(data)
    index = 0

    for i in range(len(diffs)):
        if np.all(np.abs(diffs[i : i + 500]) < threshold):
            index = i

    # Now check when the stable period ends, i.e., diff is positive > threshold
    for i in range(index, len(diffs)):
        if np.all(diffs[i : i + 500] > hysterisis*threshold):
            return (index, i)

    return False


def find_start(data, target=100):
    """
    Find the first index where data falls below target.
    Returns the index (int) or None if not found.
    """
    arr = np.asarray(data)
    idx = np.where(arr < target)[0]
    return int(idx[0]) if idx.size > 0 else None


def check_r(data, target=100, threshold=0.1):
    """
    Identifies the index in data where the value deviates from target sufficiently,
    after first approaching target from a large value.
    
    Args:
        data: Array-like data to analyze
        target: Target value to approach and deviate from
        deviation_threshold: Threshold for what constitutes sufficient deviation
        
    Returns:
        int: Index where deviation occurs, or None if pattern not found
    """
    arr = np.asarray(data)
    
    # Find when data first approaches target from above
    approach_idx = find_start(arr, target)
    if approach_idx is None:
        approach_idx = 0  # Never approached target, start from beginning
    
    # Look for deviation after approaching the target
    for i in range(approach_idx, len(arr)):
        if abs(arr[i] - target) > threshold * target:
            return i
    
    return None


def command_generation(params: dict, executable_path: str) -> str:
    """Generate data using the main.cpp file

    Args:
        params (dict): Dictionary of parameters to be fed into the main.cpp file
        executable_path (str): Path to the main.cpp executable
        errtol (float): Error tolerance
        rstep (float): Step size for robs
        pstep (float): Step size for pobs
        progress_check (int): Number of iterations to check progress, zero means no progress check

    Returns:
        str: c++ command to be executed
    """

    # Retrieve parameters to be passed as argv to the main.cpp file
    spin = params["spin"]
    charge = params["charge"]
    lam = params["lam"]
    x = params["x"]
    total_steps = params["total_steps"]
    scale = params["scale"]
    eigenswitch = params["eigenswitch"]

    # Create C++ command to be executed
    command = [
        executable_path,
        str(spin),
        str(charge),
        str(lam),
        str(x),
        str(total_steps),
        str(scale),
        str(eigenswitch),
    ]
    command = " ".join(command)

    return command


def load_geodesic(filename, cutoff=None):
    """
    Load geodesic data from C++ output file.
    """
    # Load geodesic data
    names = ['t','r', 'chi', 'phi', 'kt', 'kr', 'kchi', 'kphi', 'tau']
    data = pd.read_csv(filename, sep=' ', names=names)
    
    if cutoff:
        data = data[:cutoff]
    
    # Reverse t and tau
    # data['t'] = -data['t']
    # data['tau'] = -data['tau']
    # Remove rows with negative r
    data = data[data['r'] > 0]

    # theta = acos(chi)
    data['theta'] = np.arccos(data['chi'])
    data['ktheta'] = -data['kchi']/np.sqrt(1 - data['chi']**2)
    # Compute Cartesian coordinates
    data['x'] = data['r']*np.sin(data['theta'])*np.cos(data['phi'])
    data['y'] = data['r']*np.sin(data['theta'])*np.sin(data['phi'])
    data['z'] = data['r']*np.cos(data['theta'])
    
    return data


def load_tensor_mathematica(filename):
    """
    Load tensor data from a Mathematica output file.
    """
    # Load tidal tensor data
    names = ["eigenvalues", "eigenvectors"]
    tidal_data = pd.read_csv(filename, sep=",", names=names)

    # Convert "{A, B, C ,D}" to Python tuple
    def formater(string):
        string = string.replace("{", "[").replace("}", "]").replace("*^", "e")
        return string

    tidal_data["eigenvalues"] = tidal_data["eigenvalues"].apply(formater)
    tidal_data["eigenvectors"] = tidal_data["eigenvectors"].apply(formater)

    # Convert to float
    tidal_data["eigenvalues"] = tidal_data["eigenvalues"].apply(eval)
    tidal_data["eigenvectors"] = tidal_data["eigenvectors"].apply(eval)

    # Expand eigenvalues and eigenvectors into separate columns
    tidal_data["value_0"] = tidal_data["eigenvalues"].apply(lambda x: x[0])
    tidal_data["value_1"] = tidal_data["eigenvalues"].apply(lambda x: x[1])
    tidal_data["value_2"] = tidal_data["eigenvalues"].apply(lambda x: x[2])
    tidal_data["value_3"] = tidal_data["eigenvalues"].apply(lambda x: x[3])

    tidal_data["vector_0"] = tidal_data["eigenvectors"].apply(
        lambda x: x[0]
    )  # r vector
    tidal_data["vector_1"] = tidal_data["eigenvectors"].apply(
        lambda x: x[1]
    )  # theta vector
    tidal_data["vector_2"] = tidal_data["eigenvectors"].apply(
        lambda x: x[2]
    )  # phi vector
    tidal_data["vector_3"] = tidal_data["eigenvectors"].apply(
        lambda x: x[3]
    )  # t vector

    # For vector 0, when second component is negative, negate the vector
    def negate_vector(vector, index):
        if vector[index] < 0:
            vector = -np.array(vector)
        return vector

    tidal_data["vector_0"] = tidal_data["vector_0"].apply(negate_vector, index=1)
    tidal_data["vector_1"] = tidal_data["vector_1"].apply(negate_vector, index=2)
    tidal_data["vector_2"] = tidal_data["vector_2"].apply(negate_vector, index=3)
    tidal_data["vector_3"] = tidal_data["vector_3"].apply(negate_vector, index=0)

    return tidal_data


def load_tensor_cpp(filename, cutoff=None):
    """
    Load tensor data from a C++ output file.
    """
    # Load tidal tensor data
    names = [
        "value_0",
        "value_1",
        "value_2",
        "value_3",
        "vector_0",
        "vector_1",
        "vector_2",
        "vector_3",
    ]
    tidal_data = pd.read_csv(filename, sep=" ", names=names)
    
    if cutoff:
        tidal_data = tidal_data[:cutoff]

    # Group eigenvalues into a single column
    tidal_data["eigenvalues"] = tidal_data[["value_0", "value_1", "value_2", "value_3"]].values.tolist()

    # Rank eigenvalues from +ve to -ve
    tidal_data["eigenvalues"] = tidal_data["eigenvalues"].apply(lambda x: sorted(x, reverse=True))

    # Put eigenvalues back into separate columns
    tidal_data["value_0"] = tidal_data["eigenvalues"].apply(lambda x: x[0])
    tidal_data["value_1"] = tidal_data["eigenvalues"].apply(lambda x: x[1])
    tidal_data["value_2"] = tidal_data["eigenvalues"].apply(lambda x: x[2])
    tidal_data["value_3"] = tidal_data["eigenvalues"].apply(lambda x: x[3])
    
    # For vector rows, convert all substring "nan" to "0"
    def nan_to_zero(string):
        string = string.replace("nan", "0")
        return string
    for col in ["vector_0", "vector_1", "vector_2", "vector_3"]:
        tidal_data[col] = tidal_data[col].apply(nan_to_zero)

    # Convert vectors to Python tuples
    for col in ["vector_0", "vector_1", "vector_2", "vector_3"]:
        tidal_data[col] = tidal_data[col].apply(eval)
    
    # For vector 0, when second component is negative, negate the vector
    def negate_vector(vector, index):
        if vector[index] < 0:
            vector = -np.array(vector)
        return vector
    
    tidal_data["vector_0"] = tidal_data["vector_0"].apply(negate_vector, index=0)
    tidal_data["vector_1"] = tidal_data["vector_1"].apply(negate_vector, index=3)
    tidal_data["vector_2"] = tidal_data["vector_2"].apply(negate_vector, index=2)
    tidal_data["vector_3"] = tidal_data["vector_3"].apply(negate_vector, index=1)

    return tidal_data

def give_filename(params: dict, tidal=False) -> str:
    """
    Generate filename for the output file.
    """
    spin = params["spin"]
    charge = params["charge"]
    lam = params["lam"]
    x = params["x"]
    scale = params["scale"]
    eigenswitch = params["eigenswitch"]
    prograde = params["prograde"]

    if tidal:
        if prograde:
            filename = "data/eigensystem_a_{:.2f}_Q_{:.2f}_lambda_{:.2f}_pro.dat".format(spin, charge, lam)
        else:
            filename = "data/eigensystem_a_{:.2f}_Q_{:.2f}_lambda_{:.2f}_ret.dat".format(spin, charge, lam)
    else:
        if prograde:
            filename = "data/trace_a_{:.2f}_Q_{:.2f}_lambda_{:.2f}_pro.dat".format(spin, charge, lam)
        else:
            filename = "data/trace_a_{:.2f}_Q_{:.2f}_lambda_{:.2f}_ret.dat".format(spin, charge, lam)

    return filename

def give_params(spin, charge, lam, x, total_steps, scale, eigenswitch=1, tidal_reference=None):
    """
    Generate a dictionary of parameters.
    """
    lz = (-spin * spin + x * (x - charge * charge) - np.sqrt(x) * (x * x - 2 * x + spin * spin + charge * charge)) / (spin * (x - 1))

    if lz < 0:
        prograde = False
    else:
        prograde = True

    if eigenswitch == True:
        eigenswitch = 1
    elif eigenswitch == False:
        eigenswitch = 0

    params = {
        "spin": spin,
        "charge": charge,
        "lam": lam,
        "x": x,
        "total_steps": total_steps,
        "scale": scale,
        "eigenswitch": eigenswitch,
        "prograde": prograde,
        "lz": lz,
        "tidal_reference": tidal_reference,
    }
    
    return params