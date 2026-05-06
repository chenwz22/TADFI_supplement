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
    # The example SLiM file has been configured such that all the
    # output we want is printed on lines that start with "OUT:"
    # so we'll discard all other output lines."
    output = ""
    lines = slim_string.split('\n')
    lines = lines[24:] #去掉前面callback的
    output += lines[0][12:] + "," + lines[1][9:] + ","+ lines[2][10:] + ","+ lines[3][7:] + "," 
    # Get the specified radius / fraction / migration / growth 
    output += lines[4][9:] + "," + lines[5][7:] + ","+ lines[6][11:] + "," + lines[7][8:] + ","
    # Get the specified germline / embryo / conversion / fitness
    output += lines[8][9:] + "," + lines[9][10:] + ","
    # Get the specified dominant drive? recessive r2 drive?
    
    data = [1]
    chasing = False
    suppressed = 0
    pop_persistance = 0

    for i in range(0,len(lines)):
        if lines[i].startswith("Rates"):
            
            # When drives are newly added to the population, the population could fluctuate.
            # We mark the 5th generation's rate_dr as the baseline, and if the population rise by
            # 5% ever after, we regard it as successful invasion.
            rate_dr = lines[i].split(" ")[1][0:]
            if lines[i-2]=="generation 5" or lines[i-4]=="generation 5":   #15改成了5，我的输出里drop代=0
                data[0] = rate_dr #设定的baseline
            elif (float(rate_dr) - float(data[0])) > 0.05 and output.count(",") == 8:
                output += "1,"
        
        if lines[i].startswith("SUPPRESSED:"):
            suppressed = 1
            output = check_invasiveness(output)
            output += "1,0,0,0,0,"
            output += lines[i].split(":: ")[1] #generation
            output += ",0,0,0,0,"
            break
        elif lines[i].startswith("POP_PERSISTS:"):
            pop_persistance = 1
            output = check_invasiveness(output)
            line = lines[i-3].split(" ")[1:]   #rates
            output += "0,1,0,0,0,"
            output += lines[i].split(":: ")[1] #generation
            output += ","
            for data in line:
                output += data
                output += ","
            break

        elif lines[i].startswith("POTENTIAL_CHASE:"):
            chasing = True

        elif lines[i].startswith("LONG_TERM_CHASE:"):
            output = check_invasiveness(output)
            line = lines[i-3].split(" ")[1:]
            output += "0,0,0,0,1,"
            output += lines[i].split(":: ")[1] #generation
            output += ","
            for data in line:
                output += data
                output += ","
            break
        
    if (chasing):
        ls = check_chasing(lines)
        if (len(ls)!= 0):
            if (suppressed == 1):
                output = output.replace("1,0,0,0,0","0,0,1,0,0")
            elif (pop_persistance == 1):
                output = output.replace("0,1,0,0,0","0,0,0,1,0")
            output += str(ls["gen_chase_started"])
            output += ","
            output += str("{:.0f}".format(ls["avg_pop_during_chase"]))
            output += ","
            output += str(ls["duration_of_chasing"])
            output += ","
            output += str("{:.0f}".format(ls["avg_female_fertile"]))
        else:
            output += "0,0,0,0"
    else:
        output += "0,0,0,0"

    return output

def check_chasing(line_split):
    """
    Returns: a dictionary containing the generation when chasing starts,
    the average population during chase, the number of fertile females during chase,
    and the duration of chasing if the current simulation is a chasing condition;
    an empty dictionary if otherwise.

    This would be a chasing condition only if we find a wt allele minimum and a green's
    coefficient maximum.

    Parameter line_split: the splitted version of the entire output of a run of SLiM.
    """
    CAPACITY = 100000             #要根据capacity改
    wt_min = False
    eq_check = 0.8*2*CAPACITY
    gen_chase_started = 10000
    avg_pop_during_chase = 1000000
    duration_of_chasing = 50000
    wt = []
    gen = []
    pops = []
    gcs = []
    overall_gcs = []
    number_of_fertile_females = []
    output = {}
    for line in line_split:
        if line.startswith("FERTILE_FEMALES::"):
            spaced_line = line.split()
            number_of_fertile_females.append(int(spaced_line[1]))
        if line.startswith("WT_ALLELES::"):
            spaced_line = line.split()
            wt_alleles = int(spaced_line[1])
            this_gen = int(spaced_line[2])
            this_popsize = int(spaced_line[3])
            this_gc = float(spaced_line[5]) #gc space here
            this_overall_gc = float(spaced_line[7])
            wt.append(wt_alleles)
            gen.append(this_gen) # Creates a generations list
            pops.append(this_popsize)
            gcs.append(this_gc)
            overall_gcs.append(this_overall_gc)

    # Determine if there was a wt allele minimum
    last_gen = len(gen) - 1 # index of the last generation
    for i in range(len(wt)):
        #check 1: have at least 5 generations occurred since tracking began
        #or was this at least 5 generations prior to the end of the simulation?
        if (i > 4) and (i < (last_gen - 4)):
            this_count = wt[i]
            #check 2: is the wt count less than 80% of its eq value?
            if (this_count < eq_check):
                last_count = wt[i-1]
                next_count = wt[i+1]
                #check 3: was the last generation's wt allele count higher
                #and was the next generation's wt allele count higher?
                if (last_count > this_count) and (next_count > this_count):
                    prior_three = wt[(i-3):i]
                    next_three = wt[(i+1):(i+4)]
                    prior_avg = np.average(prior_three)
                    next_avg = np.average(next_three)
                    #check 4: was the average of the last 3 generations' wt allele
                    #counts higher and was the average of the next 3 generations'
                    #wt allele counts higher?
                    if (prior_avg > this_count) and (next_avg > this_count):
                        wt_min = True #found a minimum
                        gen_wt_min = gen[i]
                        break

    # If we've found a wt allele minimum, now check for a gc maximum
    if wt_min:
        for i in range(len(gcs)):
            #check 1: have at least 5 generations occurred since tracking began
            #or was this at least 5 generations prior to the end of the simulation?
            if (i > 4) and (i < (last_gen - 4)): #need 5 gens on each side
                this_gc_count = gcs[i]
                last_gc_count = gcs[i-1]
                next_gc_count = gcs[i+1]
                #check 2: was the last generation's green's coefficient lower
                #and was the next generation's green's coefficient count lower?
                if (last_gc_count < this_gc_count) and (next_gc_count < this_gc_count):
                    prior_three = gcs[(i-3):i]
                    next_three = gcs[(i+1):(i+4)]
                    prior_avg = np.average(prior_three)
                    next_avg = np.average(next_three)
                    #check 3: was the average of the last 3 generations' green's
                    #coefficients lower and was the average of the next 3 generations'
                    #green's coefficients lower?
                    if (prior_avg < this_gc_count) and (next_avg < this_gc_count):
                        #found both a wt_min and gc_max
                        gen_gc_max = gen[i]
                        gen_chase_started = min(gen_gc_max, gen_wt_min)
                        pos = gen.index(gen_chase_started)

                        #summary stats of chase:
                        popsizes_of_interest = pops[(pos+3):-2]
                        female_fertile_of_interest = number_of_fertile_females[(pos+3):-2]
                        avg_pop_during_chase = np.average(popsizes_of_interest)
                        #var_pop_during_chase = np.var(popsizes_of_interest)
                        avg_female_fertile = np.average(female_fertile_of_interest)
                        #overall_gcs_of_interest = overall_gcs[(pos+3):-2]
                        #overall_gc_average = np.average(overall_gcs_of_interest)
                        #overall_gc_variance = np.var(overall_gcs_of_interest)
                        duration_of_chasing = gen[-1] - gen_chase_started # index of the last generation - gen_chase_started
                        output["gen_chase_started"] = gen_chase_started
                        output["avg_pop_during_chase"] = avg_pop_during_chase
                        output["avg_female_fertile"] = avg_female_fertile
                        output["duration_of_chasing"] = duration_of_chasing
                        break
    #print(output)
    return output


def check_invasiveness(outstr):
    """
    Returns: the partially desired output.

    Add zero to output if the drive did not invade the population.
    By successful invasion, we meant a 5% frequency increase in drive rate. And
    if the drive had already invaded the populaiton, we would do nothing.

    Parameter outstr: a string containing the partially desired output
    we have generated by now.
    """
    if outstr.count(",") == 7:
        outstr += "0,"
    return outstr



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
    
    parser.add_argument('-src', '--source', default="mega_spatial_2D_homing_0703.slim", type=str,
                        help=r"SLiM file to be run. Default 'mega_spatial_2D_homing_0606.slim'")    
    
    parser.add_argument('-header', '--print_header', action='store_true', default=False,
                        help='If this is set, python prints a header for a csv file.')

    # release parameters
    parser.add_argument('-radius', '--DROP_RADIUS', default=0.2, type=float,
                        help='The drop radius of the drive. Default 0.5.')

    parser.add_argument('-introduce_ratio', '--INTRODUCE_RATIO', default=0.1, type=float,
                        help='The drop size of daisy drive . Default 0.7.')
  
    # drive parameters
    parser.add_argument('-conversion', '--DRIVE_CONVERSION', default=0.9, type=float,
                       help='Conversion of drive . Default 1.0.')
    
    parser.add_argument('-embryo_cut_rate', '--EMBRYO_RESISTANCE_CUT_RATE_F_INPUT', default=0.05, type=float,
                       help='The embryo_resistance_rate of drive . Default 1.0.')
    
    parser.add_argument('-germline_cut_rate', '--WHOLE_GERMLINE_RESISTANCE_CUT_RATE', default=0.0, type=float,
                        help='The germline_resistance_rate of drive . Default 1.0')
    
    parser.add_argument('-Dwt_female_fitness', '--SOMATIC_FITNESS_MULTIPLIER_F', default=0.8, type=float,
                       help='The embryo_resistance_rate of drive . Default 1.0.')
    
    
    # ecology parameters
    parser.add_argument('-migration', '--AVERAGE_DISTANCE', default=0.05, type=float,
                        help='The movement speed of individuals.')
    
    parser.add_argument('-growth', '--GROWTH_AT_ZERO_DENSITY', default=6.0, type=float,
                        help='The benefits that wt individuals have over the drive individuals.')

    
    # The all caps names of the following arguments must exactly match
    # the names of the constants we want to define in SLiM.


#parameters that can vary in different simulations
    #parser.add_argument('-homing', '--HOMING_SUCCESS_RATE', default=1.0, type=float,
     #                   help='The drive homing rate. Default 100 percent.')
    #parser.add_argument('-res', '--RESISTANCE_FORMATION_RATE', default=0.0, type=float,
     #                   help='The resistance formation rate. Default 0 percent.')
    #parser.add_argument('-suppression', '--RECESSIVE_FEMALE_STERILE_SUPPRESSION', action='store_true',
     #                   default=False, help='Toggles from modification drive to suppression drive.')

##add columns

    args_dict = vars(parser.parse_args())
    

    
    if args_dict.pop("print_header", None): # 需要和SLiM的输出match,代数=70
        print("radius,fraction,migration,growth,germline,embryo,conversion,fitness,dominantR2,recessiveR2"\
              "invasiveness,suppressed_withoutChase,drive_lost_withoutChase,"\
            "suppressed_afterChase,drive_lost_afterChase,long_term_chase,"\
            "generation,rate_dr,rate_has_drive,rate_r2,ending population size,"\
            "gen_chase_started,avg_pop_during_chase,duration_of_chasing,avg_female_fertile")

    # The '-header' argument prints a header for the output. This can
    # help generate a nice CSV by adding this argument to the first SLiM run:
    # Next, assemble the command line arguments in the way we want to for SLiM:
    clargs = configure_slim_command_line(args_dict)
    
#initial dataframe: first run
    slim_result = run_slim(clargs)
    #print(slim_result)

    
    parsed_result = parse_slim(slim_result)
    print(parsed_result)


    #print(parsed_result)#
if __name__ == "__main__":
    main()



