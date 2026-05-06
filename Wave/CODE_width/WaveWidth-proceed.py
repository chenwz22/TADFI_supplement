# -*- coding: utf-8 -*-
"""
Created on Sun Mar 12 14:11:57 2023

@author: Mingzuyu Pan
Modified by Weizhe
"""

from argparse import ArgumentParser
import subprocess
import numpy as np
import pandas as pd
import csv





def parse_slim(slim_string,m): 
    """
    Parse the output of SLiM to extract whatever data we're looking for.
    If we want to do a more complex analysis on the output of the SLiM file,
    this is where we do it.
    Args:
        slim_string: the entire output of a run of SLiM.
    Return
        output: the desired output we want from the SLiM simulation.
    """
    # The example SLiM file has been configured such that all the
    # output we want is printed on lines that start with "OUT:"
    # so we'll discard all other output lines."
    n=0
    position = []
    drive_carrier_frequency = []
    Num=[]
    lines = slim_string.split('\n')
    number = ['set_number']#output wave data from slim
    B_number = ['original_number']#the number of running slim
    frequency_9_1 = ['frequency_0.9_x1']
    position_9_1 = ['position_0.9_x1']
    frequency_9_2 = ['frequency_0.9_x2']
    position_9_2 = ['position_0.9_x2']
    frequency_1_1 = ['frequency_0.1_x1']
    position_1_1 = ['position_0.1_x1']
    frequency_1_2 = ['frequency_0.1_x2']
    position_1_2 = ['position_0.1_x2']
    right = ['accurate_right']
    left = ['accurate_left']
    distance = ['distance']
    
    
    
    for line in lines:
        if (line.startswith("g") == True): #generation
            n += 1  #count how many waves 
        if (line.startswith("P") == True): #POSITION                
            pos = float(line.split(":")[1]) #slice index
            drive=line.split(":")[3].strip()
            if (drive == "NULL"): # drive frequency
                drive = 10000 # no drive in this slice
            else:
                drive = float(drive) #drive freq in this slice
            position.append(pos)
            drive_carrier_frequency.append(drive)
    Num=[j for j in range(0, n) for _ in range(50)]  #50 slices for each wave
    df=pd.DataFrame([Num,position,drive_carrier_frequency],index=['Num','position','drive_carrier_frequency'])
    df = df.T
    
    #每个run里面有多个wave，提取每个wave分别计算精确width
    for i in range(n):      
        this_wave = df[df['Num'] == i]        
        # 检查Drive_frequency, 确保当前wave有>0.9的数据和<0.1的数据，否则wave无效
        Valid1 = False
        Valid2 = False
        for x in range(50):          
            if this_wave.iloc[x, 2] > 0.9: #especially for megatoxin
                Valid1 = True
            if this_wave.iloc[x, 2] < 0.1: #especially for megatoxin
                Valid2 = True     
        # print(this_wave)
        if not Valid1 or not Valid2:#wave无效
            B_number.append(m)
            number.append(i)
            right.append("NAN")
            left.append("NAN")
            frequency_9_1.append("NAN")
            position_9_1.append("NAN")
            frequency_9_2.append("NAN")
            position_9_2.append("NAN")
            frequency_1_1.append("NAN")
            position_1_1.append("NAN")
            frequency_1_2.append("NAN")
            position_1_2.append("NAN")
            distance.append('NAN')
        else:#wave有效
            #for 0.9 accuarte position
            #找到对应this_wave的index
            #index用于指示第一个<=0.9的slice位置
            index_9 = (this_wave[(this_wave['drive_carrier_frequency']<=0.9)].index.tolist())[0]
            freq_9_1 = df.iloc[index_9,2]
            posi_9_1 = df.iloc[index_9,1]
            freq_9_2 = df.iloc[(index_9-1),2]
            posi_9_2 = df.iloc[(index_9-1),1]
            frequency_9_1.append(freq_9_1)
            position_9_1.append(posi_9_1)
            frequency_9_2.append(freq_9_2)
            position_9_2.append(posi_9_2)
            # print(index_9)
            # print(freq_9_1)
            # print(posi_9_1)
            # print(freq_9_2)
            # print(posi_9_2)
            index_1 = (this_wave[(this_wave['drive_carrier_frequency']<=0.1)].index.tolist())[0]
            freq_1_1 = df.iloc[index_1,2]
            posi_1_1 = df.iloc[index_1,1]
            freq_1_2 = df.iloc[(index_1-1),2]
            posi_1_2 = df.iloc[(index_1-1),1]
            frequency_1_1.append(freq_1_1)
            position_1_1.append(posi_1_1)
            frequency_1_2.append(freq_1_2)
            position_1_2.append(posi_1_2)
            # print(index_1)
            # print(freq_1_1)
            # print(posi_1_1)
            # print(freq_1_2)
            # print(posi_1_2)
            accurate_1 = posi_1_2 + (freq_1_2-0.1)/(freq_1_2-freq_1_1)*0.02
            accurate_1 = round(accurate_1,6)
            # print(accurate_1)
            accurate_9 = posi_9_2 + (freq_9_2-0.9)/(freq_9_2-freq_9_1)*0.02
            accurate_9 = round(accurate_9,6)
            # print(accurate_9)
            distance_now = accurate_1 - accurate_9
            distance_now = round(distance_now,6)
            if (posi_9_1 == 0.01):
                distance_now = "NAN"
            if (freq_9_1 <= freq_1_1): 
                distance_now = "NAN" 
            if (freq_9_2 == 0.0):# allowing freq_9_2==1
                distance_now = "NAN"                
            if (freq_1_2-freq_1_1 == 0):
                accurate_1 = 'NAN'
            if (freq_9_2-freq_9_1==0):
                accurate_9 = 'NAN'
            number.append(i)
            B_number.append(m) #####
            right.append(accurate_1)
            left.append(accurate_9)
            distance.append(distance_now)
    rows = pd.DataFrame(zip (B_number,number,position_9_2,frequency_9_2,position_9_1,frequency_9_1,position_1_2,frequency_1_2,position_1_1,frequency_1_1,left,right,distance))
    # print(number)
    # print(right)
    # print(left)
    # print(distance)
    return rows


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
    #print(err)
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


def main(replicates):
    """
    1. Configure using argparse.
    2. Generate the command line list to pass to subprocess through the run_slim() function.
    3. Run SLiM.
    4. Process the output of SLiM to extract the information we want.
    5. Print the results.
    """
    # Get args from arg parser:
    for x in range(replicates):
        parser = ArgumentParser()
        
        parser.add_argument('-src', '--source', default="homingdom_width.slim", type=str,
                            help=r"SLiM file to be run. Default 'megaTOXIN_Drive.slim'")
        
        parser.add_argument('-header', '--print_header', action='store_true', default=False,
                            help='If this is set, python prints a header for a csv file.')
      
        parser.add_argument('-embryo_resistance_rate', '--EMBRYO_RESISTANCE_RATE', default=0.1, type=float,
                           help='The embryo_resistance_rate of drive . Default 1.0.')
        
        parser.add_argument('-germline_resistance_rate', '--GERMLINE_RESISTANCE_RATE', default=0.1, type=float,
                            help='The germline_resistance_rate of drive . Default 1.0')
        
        parser.add_argument('-conversion', '--GERMLINE_EFFICIENCY', default=0.9, type=float,
                            help='The germline_resistance_rate of drive . Default 1.0')
        
        parser.add_argument('-fitness', '--SOMATIC_FITNESS_MUTLIPLIER_F', default=0.9, type=float,
                            help='The germline_resistance_rate of drive . Default 1.0')
        #parser.add_argument('-capacity', '--CAPACITY', default=100000, type=int,
        #                    help='The germline_resistance_rate of drive . Default 1.0')
        
        #parser.add_argument('-drop_size', '--DROP_SIZE', default=10000, type=int,
        #                    help='The germline_resistance_rate of drive . Default 1.0')
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
        # print(slim_result)
        #with open("mega_wavewidth0824.txt", "w", encoding="utf-8") as f:
        #        f.write(slim_result)
        parsed_result = parse_slim(slim_result,x)    
        # print(parsed_result)
        parsed_result.to_csv("homingdom_wavewidth0825.csv",encoding = "gbk",index=False,header=False,mode='a', )
        
    
if __name__ == "__main__":
    main(20)



