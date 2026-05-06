import subprocess
import numpy as np
from argparse import ArgumentParser
import re

def parse_slim(slim_string):
    """
    Parse the output of SLiM to extract whatever data we're looking for.
    """
    lines = slim_string.split('\n')
    lines = lines[24:] # 去掉前面callback的
    
    output = ""
    
    # 修改参数映射表以匹配SLiM的实际输出
    parameter_map = {
        "DROP_RADIUS:": "radius",
        "FRACTION:": "introduce_ratio",
        "Migration:": "migration",
        "R1_Occurrence:": "r1_occurrence",
        "Germline:": "germline",
        "Embryo:": "embryo",
        "Conversion:": "conversion",
        "Fitness:": "fitness",
        "Dominant:": "dominant",
        "Recessive:": "recessive"
    }
    
    # 初始化参数字典
    params = {value: "unknown" for value in parameter_map.values()}
    
    # 扫描前20行查找参数
    for i in range(min(20, len(lines))):
        line = lines[i].strip()
        for prefix, param_name in parameter_map.items():
            if line.startswith(prefix):
                # 提取参数值
                value = line[len(prefix):].strip()
                params[param_name] = value
    
    # 特殊处理：将r1_occurrence映射到growth
    if "r1_occurrence" in params and params["r1_occurrence"] != "unknown":
        params["growth"] = params["r1_occurrence"]
    
    # 按顺序构建输出字符串
    param_order = ["radius", "introduce_ratio", "migration", "growth", 
                  "germline", "embryo", "conversion", "fitness", 
                  "dominant", "recessive"]
    
    for param in param_order:
        output += params[param] + ","
    
    # 初始化变量
    data = [1]
    chasing = False
    suppressed = 0
    pop_persistance = 0
    pop_resistence = 0
    is_long_term_chase = False
    long_term_avg_female_fertile = None
    is_pop_resistence = False
    avg_female_fertile_resistence = None
    number_of_fertile_females = []
    fertile_females_all = []
    for line in lines:
        if line.startswith("FERTILE_FEMALES::"):
            parts = line.split()
            if len(parts) > 1:
                try:
                    fertile_females_all.append(int(parts[1]))
                except ValueError:
                    pass
    final_female_fertile = str(fertile_females_all[-1])

    # 先收集所有的 FERTILE_FEMALES 数据
    for line in lines:
        if line.startswith("FERTILE_FEMALES::"):
            try:
                number_of_fertile_females.append(int(line.split()[1]))
            except (IndexError, ValueError):
                pass
    
    # 继续处理其他数据
    for i in range(0, len(lines)):
        if i < len(lines) and lines[i].startswith("Rates"):
            try:
                rate_parts = lines[i].split(" ")
                if len(rate_parts) > 1:
                    rate_dr = rate_parts[1]
                    
                    if (i-2 < len(lines) and lines[i-2]=="generation 5") or (i-4 < len(lines) and lines[i-4]=="generation 5"):
                        data[0] = rate_dr # 设定的baseline
                    elif float(rate_dr) - float(data[0]) > 0.05 and output.count(",") == 10:
                        output += "1,"
            except Exception:
                pass
        
        if lines[i].startswith("SUPPRESSED:"):
            suppressed = 1
            output = check_invasiveness(output)
            output += "1,0,0,0,0,0,0,"
            try:
                output += lines[i].split(":: ")[1] # generation
            except IndexError:
                output += "unknown,"
            output += ",0,0,0,0,"
            break
            
        elif lines[i].startswith("POP_PERSISTS:"):
            pop_persistance = 1
            output = check_invasiveness(output)
            try:
                line = lines[i-3].split(" ")[1:]   # rates
                output += "0,1,0,0,0,0,0,"
                output += lines[i].split(":: ")[1] # generation
                output += ","
                for data in line:
                    output += data + ","
            except Exception:
                output += "0,1,0,0,0,0,0,unknown,unknown,unknown,unknown,unknown,"
            break

        elif lines[i].startswith("RESISTANCE:"):
            pop_resistence = 1
            output = check_invasiveness(output)
            
            try:
                line = lines[i-3].split(" ")[1:]  # rates
                output += "0,0,1,0,0,0,0,"  # new category
                output += lines[i].split(":: ")[1]  # generation
                output += ","
                for data in line:
                    output += data + ","
            except Exception:
                output += "0,0,1,0,0,0,0,unknown,unknown,unknown,unknown,unknown,"
            
            # 计算平均fertile female
            if number_of_fertile_females:
                avg_female_fertile_resistence = int(round(np.average(number_of_fertile_females)))
            else:
                avg_female_fertile_resistence = "not detected"
            
            is_pop_resistence = True
            break

        elif lines[i].startswith("POTENTIAL_CHASE:"):
            chasing = True

        elif lines[i].startswith("LONG_TERM_CHASE:"):
            output = check_invasiveness(output)
            try:
                line = lines[i-3].split(" ")[1:]
                output += "0,0,0,0,0,0,1,"
                output += lines[i].split(":: ")[1]  # generation
                output += ","
                for data in line:
                    output += data + ","
            except Exception:
                output += "0,0,0,0,0,0,1,unknown,unknown,unknown,unknown,unknown,"
                
            fertile_females = []
            for l in lines:
                if l.startswith("FERTILE_FEMALES::"):
                    try:
                        fertile_females.append(int(l.split()[1]))
                    except (IndexError, ValueError):
                        pass
            
            if fertile_females:
                long_term_avg_female_fertile = int(round(sum(fertile_females) / len(fertile_females)))
            else:
                long_term_avg_female_fertile = "not detected"
                
            is_long_term_chase = True
            break

    if (chasing):
        ls = check_chasing(lines)
        if (len(ls) != 0):
            if (suppressed == 1):
                output = output.replace("1,0,0,0,0,0,0,", "0,0,0,1,0,0,0,")
            elif (pop_persistance == 1):
                output = output.replace("0,1,0,0,0,0,0,", "0,0,0,0,1,0,0,")
            elif (pop_resistence == 1):
                output = output.replace("0,0,1,0,0,0,0,", "0,0,0,0,0,1,0,")
            output += str(ls["gen_chase_started"]) + ","
            output += str("{:.0f}".format(ls["avg_pop_during_chase"])) + ","
            output += str(ls["duration_of_chasing"]) + ","
            output += str("{:.0f}".format(ls["avg_female_fertile"]))
        else:
            if is_long_term_chase and long_term_avg_female_fertile is not None:
                output += "0,0,0,0"
            elif is_pop_resistence and avg_female_fertile_resistence is not None: 
                output += "0,0,0,0"
            else:
                output += "0,0,0,0"
    else:
        if is_long_term_chase and long_term_avg_female_fertile is not None:
            output += "0,0,0,0"
        elif is_pop_resistence and avg_female_fertile_resistence is not None: 
            output += "0,0,0,0"
        else:
            output += "0,0,0,0"
    output += "," + final_female_fertile
    return output

def check_chasing(line_split):
    """
    Returns: a dictionary containing the generation when chasing starts,
    the average population during chase, the number of fertile females during chase,
    and the duration of chasing if the current simulation is a chasing condition;
    an empty dictionary if otherwise.
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
            try:
                spaced_line = line.split()
                if len(spaced_line) > 1:
                    number_of_fertile_females.append(int(spaced_line[1]))
            except (IndexError, ValueError):
                pass
                
        if line.startswith("WT_ALLELES::"):
            try:
                spaced_line = line.split()
                if len(spaced_line) >= 8:
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
            except (IndexError, ValueError):
                pass

    # Determine if there was a wt allele minimum
    if len(gen) > 0:
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
                            if pos + 3 < len(pops) and len(pops) > 2:
                                popsizes_of_interest = pops[(pos+3):-2]
                                if len(number_of_fertile_females) > pos + 3 and len(number_of_fertile_females) > 2:
                                    female_fertile_of_interest = number_of_fertile_females[(pos+3):-2]
                                    avg_pop_during_chase = np.average(popsizes_of_interest)
                                    avg_female_fertile = np.average(female_fertile_of_interest)
                                    duration_of_chasing = gen[-1] - gen_chase_started
                                    output["gen_chase_started"] = gen_chase_started
                                    output["avg_pop_during_chase"] = avg_pop_during_chase
                                    output["avg_female_fertile"] = avg_female_fertile
                                    output["duration_of_chasing"] = duration_of_chasing
                                    break
    
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
    if outstr.count(",") == 10:  # 10个初始参数
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
    
    parser.add_argument('-src', '--source', default="mega_spatial_2D_homing_0703_dom_LJY.slim", type=str,
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
    
    parser.add_argument('-embryo_cut_rate', '--EMBRYO_RESISTANCE_CUT_RATE_F_INPUT', default=0.03, type=float,
                       help='The embryo_resistance_rate of drive . Default 1.0.')
    
    parser.add_argument('-germline_cut_rate', '--WHOLE_GERMLINE_RESISTANCE_CUT_RATE', default=0.1, type=float,
                        help='The germline_resistance_rate of drive . Default 1.0')
    
    parser.add_argument('-Dwt_female_fitness', '--SOMATIC_FITNESS_MULTIPLIER_F', default=0.8, type=float,
                       help='The embryo_resistance_rate of drive . Default 1.0.')

    parser.add_argument('-r1_occurrence', '--R1_OCCURRENCE_RATE', default=0.14, type=float,
                       help='The embryo_resistance_rate of drive . Default 1.0.')
    
    # ecology parameters
    parser.add_argument('-migration', '--AVERAGE_DISTANCE', default=0.05, type=float,
                        help='The movement speed of individuals.')
    
    parser.add_argument('-growth', '--GROWTH_AT_ZERO_DENSITY', default=6.0, type=float,
                        help='The benefits that wt individuals have over the drive individuals.')

    args_dict = vars(parser.parse_args())
    
    # 如果需要打印表头
    if args_dict.pop("print_header", None): # 需要和SLiM的输出match,代数=70
        print("radius,fraction,migration,r1_occurrence,germline,embryo,conversion,fitness,dominantR2,recessiveR2,"\
              "invasiveness,suppressed_withoutChase,drive_lost_withoutChase,pop_resistence_without_chase,"\
            "suppressed_afterChase,drive_lost_afterChase,pop_resistence_after_chase, long_term_chase,"\
            "generation,rate_dr,rate_has_drive,rate_r2,ending_population_size,"\
            "gen_chase_started,avg_pop_during_chase,duration_of_chasing,avg_female_fertile,final_female_fertile")

    # 配置SLiM命令行参数
    command_line_args = configure_slim_command_line(args_dict)
    
    # 运行SLiM
    slim_output = run_slim(command_line_args)
    
    # 解析SLiM输出
    parsed_output = parse_slim(slim_output)
    
    # 打印结果
    print(parsed_output)


if __name__ == "__main__":
    main()
