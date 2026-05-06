#!/usr/bin/env python3

#  Created by Sam Champer, 2020.
#  A product of the Messer Lab, http://messerlab.org/slim/

#  Sam Champer, Ben Haller and Philipp Messer, the authors of this code, hereby
#  place the code in this file into the public domain without restriction.
#  If you use this code, please credit SLiM-Extras and provide a link to
#  the SLiM-Extras repository at https://github.com/MesserLab/SLiM-Extras.
#  Thank you.

# This is an example of how to use Python as a driver for SLiM.
# Output is formated as a csv, but just printed to stdout by default.

# In order to reconfigure this file for your research project, the
# run_slim() and configure_slim_command_line() functions do not need to be modified.
# Changes you would likely want to make are to the argument parser in main(),
# in order to pass your desired variables to SLiM, and to the parse_slim()
# function, where you could do your desired operations on the output of SLiM.



#try: 10 times each parameter
#divide a list into several sub_lists


from argparse import ArgumentParser
import subprocess
import numpy as np
from math import ceil
import csv
import pandas as pd

def parse_slim(slim_string,m):
    """
    Parse the output of SLiM to extract whatever data we're looking for.
    If we want to do a more complex analysis on the output of the SLiM file,
    this is where we do it.
    Args:
        slim_string: the entire output of a run of SLiM.
    Return
        output: the desired output we want from the SLiM simulation.
    Tocsv
        output parsed result to csv
    """
    lines = slim_string.strip().splitlines()
    
    data = []
    row = {}
    for i, line in enumerate(lines):
        if "OUT::" in line:
            lines = lines[(i+1):]  # 从"OUT::"之后开始保留
            break
    for line in lines: 
        if line == "":
            continue  # skip blank lines
        if line.startswith("Gen:"):
            if row:
                data.append(row)
            row = {"Gen": line.split(":")[1].strip()}
        else:
            key, value = line.split(":",1)
            row[key.strip()] = value.strip()
    if row:
        data.append(row)
    fieldnames = sorted({key for d in data for key in d.keys()})
    filename="".join(["DOMINANT_panmictic", str(m), ".csv"])    #记得改文件名！
    
    with open(filename, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(data)
    return data



def run_slim(command_line_args):
    """
    Runs SLiM using subprocess.
    Args:
        command_line_args: list; a list of command line arguments.
    return: The entire SLiM output as a string.
    """
    slim = subprocess.Popen(command_line_args,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            universal_newlines=True)
    out, err = slim.communicate()
    # For debugging purposes:
    # std.out from the subprocess is in slim.communicate()[0]
    # std.error from the subprocess is in slim.communicate()[1]
    # Errors from the process can be printed with:
    # print(err)
    return out


def configure_slim_command_line(args_dict):
    """
    Sets up a list of command line arguments for running SLiM.
    Args:
        args_dict: a dictionary of arg parser arguments.
    Return
        clargs: A formated list of the arguments.
    """
    # We're running SLiM, so the first arg is simple:
    clargs = "slim "
    # The filename of the source file must be the last argument:
    source = args_dict.pop("source")
    # Add each argument from arg parser to the command line arguemnts for SLiM:
    for arg in args_dict:
        if isinstance(args_dict[arg], bool):
            clargs += f"-d {arg}={'T' if args_dict[arg] else 'F'} "
        else:
            clargs += f"-d {arg}={args_dict[arg]} "
    # Add the source file, and return the string split into a list.
    clargs += source
    return clargs.split()


def main(x):
    """
    1. Configure using argparse.
    2. Generate the command line list to pass to subprocess through the run_slim() function.
    3. Run SLiM.
    4. Process the output of SLiM to extract the information we want.
    5. Print the results.
    """
    # Get args from arg parser:
    
    for i in range(x):
        parser = ArgumentParser()
        parser.add_argument('-src', '--source', default="mega_panmictic_0612.slim", type=str,
                            help=r"SLiM file to be run. Default 'megaTOXIN_Drive.slim'")
        
        
        parser.add_argument('-header', '--print_header', action='store_true', default=False,
                            help='If this is set, python prints a header for a csv file.')
    
        parser.add_argument('-introduce_ratio', '--INTRODUCE_RATIO', default=0.1, type=float,
                            help='The drop size of daisy drive . Default 0.7.')
      
        parser.add_argument('-embryo_resistance_rate', '--DRIVE_HETEROZYGOTE_EMBRYO_RESISTANCE_RATE', default=0.1, type=float,
                           help='The embryo_resistance_rate of drive . Default 1.0.')
        
        parser.add_argument('-germline_resistance_rate', '--GERMLINE_RESISTANCE_RATE', default=0.2, type=float,
                            help='The germline_resistance_rate of drive . Default 1.0')
        
        parser.add_argument('-drive_conversion', '--DRIVE_CONVERSION_RATE', default=0.6, type=float,
                            help='The germline_resistance_rate of drive . Default 1.0')
        
        parser.add_argument('-Dwt_female_fitness', '--SOMATIC_FITNESS_MUTLIPLIER_F', default=0.8, type=float,
                           help='The embryo_resistance_rate of drive . Default 1.0.')
    
    
    
    # The all caps names of the following arguments must exactly match
    # the names of the constants we want to define in SLiM.



##add columns

        args_dict = vars(parser.parse_args())
    
    
        if args_dict.pop("print_header", None): # 需要和SLiM的输出match,代数=70
            print("generation,introduce_ratio,drive_conversion,germline_resistance,embryo_resistance,Dwt_F_fitness,pop_size,rate_wt,rate_dr,rate_complete_r2,success_Yes,")

    # The '-header' argument prints a header for the output. This can
    # help generate a nice CSV by adding this argument to the first SLiM run:
    # Next, assemble the command line arguments in the way we want to for SLiM:
        clargs = configure_slim_command_line(args_dict)

#initial dataframe: first run
        slim_result = run_slim(clargs)
    
        parsed_result = parse_slim(slim_result,i+1)

        #print(parsed_result)


    #print(parsed_result)#
if __name__ == "__main__":
    main(20)



