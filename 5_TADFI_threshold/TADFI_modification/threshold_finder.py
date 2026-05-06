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


from argparse import ArgumentParser
import subprocess


def parse_slim(slim_string):
    """
    Parse the output of SLiM to extract whatever data we're looking for.
    If we want to do a more complex analysis on the output of the SLiM file,
    this is where we do it.
    Args:
        slim_string: the entire output of a run of SLiM.
    Return
        output: the desired output we want from the SLiM simulation.
    """
    lines = slim_string.split('\n')
    for line in lines:
        if line.startswith("OUT:"):
            contents = line[5:].split(",")
            above_threshold = float(contents[-2])            
    return above_threshold


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
    args_dict["source"] = source
    return clargs.split()


def main():
    """
    1. Configure using argparse.
    2. Generate the command line list to pass to subprocess through the run_slim() function.
    3. Run SLiM.
    4. Process the output of SLiM to extract the information we want.
    5. Print the results.
    """
    # Get args from arg parser:
    parser = ArgumentParser()
    parser.add_argument('-src', '--source', default="megaTOXIN_Drive_TH2_modification.slim", type=str, ##########记得改回来！！！！！！！1
                        help=r"SLiM file to be run. Default 'megaTOXIN_Drive_TH2.slim'")
    parser.add_argument('-header', '--print_header', action='store_true', default=False,
                        help='If this is set, python prints a header for a csv file.')
    parser.add_argument('-embryo', '--DRIVE_HETEROZYGOTE_EMBRYO_RESISTANCE_RATE', default=0.0, type=float)
    parser.add_argument('-germline', '--GERMLINE_RESISTANCE_RATE', default=1.0, type=float)
    parser.add_argument('-DD_fitness', '--DD_FITNESS_VALUE', default=1.0, type=float)
    

    args_dict = vars(parser.parse_args())

    if args_dict.pop("print_header", None):
        # The '-header' argument prints a header for the output. This can
        # help generate a nice CSV by adding this argument to the first SLiM run:
        # Print the variable names.
        #print(','.join(f"{arg}" for arg in args_dict if arg != "source"), end=",")
        # Print a heading for the data being collected:
        #print("Minimum introduction threshold")
        print("embryo_cut_rate,germline_cut_rate,DD_fitness,threshold,capacity,")

    ####################
    ####################
    # Test introduction frequency for 5 times, if "Nearly Threshold" < 4 times or succeed < 2 times
    # output threshold=1.0, no more binary search and accurate search later.    
    args_dict["INTRODUCE_RATIO"] = 0.99   
    PRE_current_level_above_threshold = 0
    for _ in range(3):
        args_dict["CAPACITY"] = 20000  ####
        clargs = configure_slim_command_line(args_dict)
        slim_result = run_slim(clargs)
        PRE_current_level_above_threshold += parse_slim(slim_result)
        #print(PRE_current_level_above_threshold)
    if PRE_current_level_above_threshold < 1:    #降低了进入循环的要求，1/3成功就能进入          
        # NO more try
        args_dict["INTRODUCE_RATIO"] = 1.0
    else:        
        #continue binary search
        args_dict["CAPACITY"] = 60000 #change from 20000 to 60000
        args_dict["INTRODUCE_RATIO"] = 0.5
        step_size = 0.25
        # Search for the required introduction frequency with a binary search.   
        for _ in range(6):
            # Binary search will be go through 6 introduction ratios.    
            # To decide whether to go up or down, we check how many
            # simulations are above the threshold at this current introduction ratio.
            current_level_above_threshold = 0
            sims_per_level = 5
            for _ in range(sims_per_level):
                #print(args_dict["INTRODUCE_RATIO"])########
                # Run 5 simulations at the current introduction frequency.
                # If the drive usually fails, we'll increase the frequency.
                # If the drive usually succeeds, we'll decrease the frequency,
    
                # Assemble the command line arguments in the way we want to for SLiM:
                clargs = configure_slim_command_line(args_dict)
                slim_result = run_slim(clargs)
                current_level_above_threshold += parse_slim(slim_result)
            if current_level_above_threshold > sims_per_level / 2:
                # Since the SLiM code returns 0.5 when the drive is near threhold,
                # I expect that the drive is above it's threshold in this case.
                # Decrease the introduction ratio, and divide step size by two.
                args_dict["INTRODUCE_RATIO"] -= step_size
                step_size /= 2
            else:
                # Otherwise, increase the introduction rate.
                args_dict["INTRODUCE_RATIO"] += step_size
                step_size /= 2
        
        # Next, verify the threshold using a linear search and more simulation replicates.
        # Start with introduction ratio 5% lower and truncated to the nearest 100th.        
        args_dict["CAPACITY"] = 100000
        args_dict["INTRODUCE_RATIO"] = float(f"{args_dict['INTRODUCE_RATIO'] -0.05 :-0.2f}") #5% lower
        if  args_dict["INTRODUCE_RATIO"]< 0.0: #避免负数
            args_dict["INTRODUCE_RATIO"]=0.0
        while True:     
            current_level_above_threshold = 0          
            sims_per_level = 10
            for _ in range(sims_per_level):
                #print(args_dict["INTRODUCE_RATIO"])############
                clargs = configure_slim_command_line(args_dict)
                slim_result = run_slim(clargs)
                current_level_above_threshold += parse_slim(slim_result)
            #print(current_level_above_threshold)
            if current_level_above_threshold < sims_per_level / 2 and args_dict["INTRODUCE_RATIO"] < 0.99:
                # Increase by 1 percent.
                args_dict["INTRODUCE_RATIO"] = float(f"{args_dict['INTRODUCE_RATIO'] +0.01 :-0.2f}")   #######to avoid 0.37000000000000005
            elif current_level_above_threshold < sims_per_level / 2 and args_dict["INTRODUCE_RATIO"] == 0.99:
                # assume never succeed in any introduction_ratio<=0.99, print 1.0
                args_dict["INTRODUCE_RATIO"] = 1.0
                break
            else:
                # This run was at the threshold, so we're done.
                break

    # This next line will print out all of the paramters for this simulation,
    # along with the introduction ratio that was determined to be the minimum threshold.
    print(','.join(f"{args_dict[arg]}" for arg in args_dict if arg != "source"), end=",\n")
    

    ####################
    ####################

if __name__ == "__main__":
    main()
