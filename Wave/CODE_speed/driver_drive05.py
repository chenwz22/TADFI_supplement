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

#在Pan的代码里面有两种情况，一种是用drive>0.5进行判断，另一种是用wt<0.8进行判断
#这两种情况的slim输出和python数据处理是需要配套用的
#我的megatoxin属于wt<0.8的情况；slim输出的是wt—freq，python需要处理的也是针对wt-freq的
#slice-helper属于记录变量，相当于文献里公式的G0或F0
#我写了两份代码，本代码是给drive>0.5的

from argparse import ArgumentParser
import subprocess
import pandas as pd

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
   
    SLICE3 =['SLICE3']
    SLICE8 = ['SLICE8']

    Before_Start_Generation = ['Before_Start_Generation']
    Before_Start = ['Before_Start']
    Start_Generation = ['Start_Generation']
    Start = ['Start']
    Actual_Start_Generation = ['Actual_Start_Generation']

    Before_Stop_Generation = ['Before_Stop_Generation']
    Before_Stop = ['Before_Stop']
    Stop_Generation = ['Stop_Generation']
    Stop = ['Stop']
    Actual_Stop_Generation = ['Actual_Stop_Generation']
    Actual_Timed_Gens = ['Actual_Timed_Gens']

    Embryo = ['EMBRYO']
    Germline = ['GERMLINE']
    Number = ['Number']  
    Timed_Gens_Tester = ['Timed_Gens_Tester']  
    #Calculate = 'NULL'
    MAX = 500
        
    
        
    lines = slim_string.split('\n')
    #print(lines)
    slice3_helper = 0
    slice8_helper = 0
    actual_start = 'NULL'
    actual_stop = 'NULL'
    for line in lines:
        if line.startswith("generation"):
            gen = line.split(':')[1].strip()                
        if line.startswith("SLICE3"):
            #print(line)
            SLICE3_1 = line.split(':')[1].strip()                
            if(SLICE3_1 == 'N/A'):
                SLICE3_1 = '10'
            SLICE3_1 = float(SLICE3_1)
            SLICE3.append(SLICE3_1)
            
            if (SLICE3_1 < 0.5) & (slice3_helper!= "NULL"):
                slice3_helper = SLICE3_1            
            if (SLICE3_1 >= 0.5) &(slice3_helper!= "NULL") :            
                Before_Start.append(slice3_helper) 
                Before_Start_Generation.append(float(gen)-1) 
                Start.append(SLICE3_1)
                Start_Generation.append(float(gen))
                actual_start = float(gen)-1 + (0.5 - float(slice3_helper))/(SLICE3_1 - float(slice3_helper))
                Actual_Start_Generation.append(actual_start)
                slice3_helper = "NULL" 
                
                
        if line.startswith("SLICE8"):
            SLICE8_0 = line.split(':')[1]
            SLICE8_1 = SLICE8_0.split('\n')[0]
            float_SLICE8_1 = float(SLICE8_1)
            SLICE8.append(SLICE8_1)
            if (float_SLICE8_1 < 0.5) & (slice8_helper!= "NULL"):
                slice8_helper = SLICE8_1
            if (slice8_helper!= "NULL") & (float_SLICE8_1 >= 0.5):
                Before_Stop.append(slice8_helper)
                Before_Stop_Generation.append(float(gen)-1)
                Stop.append(float_SLICE8_1)
                Stop_Generation.append(float(gen))
                actual_stop = float(gen)-1 + (0.5 - float(slice8_helper))/(-float(slice8_helper) + float_SLICE8_1)
                Actual_Stop_Generation.append(actual_stop)
                slice8_helper = "NULL" 
                
        if line.startswith("EMBRYO"):
            Embryo_1 = line.split(':')[1].strip()
            Embryo.append(Embryo_1)
        if line.startswith("GERMLINE"):
            Germline_1 = line.split(':')[1].strip()
            Germline.append(Germline_1)
            Number.append(1) ##如果需要对i个输出求平均，此处改成i，遍历文件夹所有输出文件
            
        if line.startswith("TIMED"):
            Timed_Gens_Tester_1 = line.split(':')[1].strip()
            Timed_Gens_Tester.append(Timed_Gens_Tester_1)
      
            if((actual_start != 'NULL') & (actual_stop != 'NULL')): #正常来说进入这个循环
                    actual_timed = actual_stop-actual_start
                    Actual_Timed_Gens.append(actual_timed)                
        
            if((actual_stop == 'NULL') & (actual_start == 'NULL')): 
                    Before_Start_Generation.append('0')
                    Before_Start.append('0')
                    Start_Generation.append('0')
                    Start.append('0')
                    Before_Stop_Generation.append('0')
                    Before_Stop.append('0')
                    Stop_Generation.append('0')
                    Stop.append('0')
                    Actual_Start_Generation.append(MAX)
                    Actual_Stop_Generation.append(MAX)
                    Actual_Timed_Gens.append(MAX)

            if((actual_stop == 'NULL') & (actual_start != 'NULL')):
                    Before_Stop_Generation.append('0')
                    Before_Stop.append('0')
                    Stop_Generation.append('0')
                    Stop.append('0')
                    Actual_Stop_Generation.append(MAX)
                    Actual_Timed_Gens.append(MAX)

            if((actual_stop != 'NULL') & (actual_start == 'NULL')):
                    Before_Start_Generation.append('0')
                    Before_Start.append('0')
                    Start_Generation.append('0')
                    Start.append('0')
                    Actual_Start_Generation.append(MAX)
                    Actual_Timed_Gens.append(MAX)
    rows = zip(Number,Germline,Embryo,Before_Start_Generation,Before_Start,Start_Generation,Start,Actual_Start_Generation,Before_Stop_Generation,Before_Stop,Stop_Generation,Stop,Actual_Stop_Generation,Timed_Gens_Tester,Actual_Timed_Gens)
    df = pd.DataFrame(rows)
    #print(Before_Start_Generation)
    return df


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
    parser.add_argument('-src', '--source', default="homingdomDD_speed.slim", type=str, ##########记得改回来！！！！！！！1
                        help=r"SLiM file to be run. Default 'homingdom_speed.slim'")
    parser.add_argument('-header', '--print_header', action='store_true', default=False,
                        help='If this is set, python prints a header for a csv file.')
    parser.add_argument('-embryo', '--EMBRYO_RESISTANCE_RATE', default=0.1, type=float)
    parser.add_argument('-germline', '--GERMLINE_RESISTANCE_RATE', default=0.1, type=float)
    parser.add_argument('-fitness', '--SOMATIC_FITNESS_MUTLIPLIER_F', default=0.9, type=float)
    parser.add_argument('-homing', '--GERMLINE_EFFICIENCY', default=0.9, type=float)

    args_dict = vars(parser.parse_args())

    if args_dict.pop("print_header", None):
        # The '-header' argument prints a header for the output. This can
        # help generate a nice CSV by adding this argument to the first SLiM run:
        # Print the variable names.
        print("embryo_cut_rate,germline_cut_rate,Dwt_female_fitness,homing,threshold,capacity,")

        print(','.join(f"{args_dict[arg]}" for arg in args_dict if arg != "source"), end=",\n")

     # Next, assemble the command line arguments in the way we want to for SLiM:
    clargs = configure_slim_command_line(args_dict)
    
     # Run the file with the desired arguments.
    slim_result = run_slim(clargs)
       
     # Parse and analyze the result.
    parsed_result = parse_slim(slim_result)
    #print(parsed_result)
    parsed_result.to_csv("homingdomDD_speed0901.csv",encoding = "gbk",index=False,header=False,mode='a',)

if __name__ == "__main__":
    main()
