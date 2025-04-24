#!/usr/bin/env python
import os,sys
import time
import logging
import argparse
from itertools import islice
from multiprocessing import Pool,cpu_count 
from multiprocessing.dummy import Pool as ThreadPool


class HelpFormatter(argparse.RawDescriptionHelpFormatter,argparse.ArgumentDefaultsHelpFormatter):
    pass

def main():
    args      = get_args()
    inputfile = check_input(args.INPUT)
    workdir   = os.path.dirname(inputfile)
    logfile   = inputfile+".log"
    rundir    = inputfile+".run"

    if args.f :
        print("#Clean the run dir and log file ")
        os.system( "rm -rf %s.* " %(inputfile) )
    else:
        check_done( "%s.done"%(os.path.abspath(args.INPUT)) ,"break")

    if args.m > cpu_count()-2 :
        args.m = cpu_count()-2



    job_list=creat_sub_cmd(inputfile,args.l,args.n)

    if args.t == "local":
        pool = ThreadPool(args.m) 
        pool.map(local_run,job_list)
        pool.close()
        pool.join()
    elif args.t == "pbs":
        cluster_run(job_list,args.ncpu,args.queue,"pbs")
    elif args.t == "sge":
        cluster_run(job_list,args.ncpu,args.queue,"sge")
    elif args.t == "slurm":
        cluster_run(job_list,args.ncpu,args.queue,"slurm")
    else:
        print("# Error! This program only support SGE PBS SLURM and Local to run commond")
        exit(1)


    check_all(inputfile,job_list)


def runlocal(inputfile,lines=1,threads=8):

    dones=check_done( "%s.done"%(os.path.abspath(inputfile)) ,"continue")
    if dones == 1:
        print(f"# Attention cmd sikpped! due to {os.path.abspath(inputfile)} exisit")
        return 1

    logfile   = inputfile+".log"
    logging.basicConfig(level=logging.DEBUG,  
            format='%(asctime)s  %(levelname)s %(message)s',  
            datefmt='%a, %d %b %Y %H:%M:%S',  
            filename=logfile,  
            filemode='a')      

    job_list=creat_sub_cmd(inputfile,lines,"work")

    pool = ThreadPool(threads) 
    pool.map(local_run,job_list)
    pool.close()
    pool.join()
    check_all(inputfile,job_list)


def check_all(inputfile,job_list):
    finish = 0
    undone = ""
    for job in job_list :
        if os.path.exists(f"{job}.done"):
            finish += 1 
        else:
            undone += f"{job}\n"
    
    if finish < len(job_list):
        print(f"Failed! only {finish} out of {len(job_list)} is done!")
        print(f"Please check your commond: {inputfile}")
        logging.error( "Failed! only %s out of %s is done!" %(finish,len(job_list)) )
        logging.error( "Please check your undone list: %s.undone"% (inputfile)  )
        with open(f"{inputfile}.undone","w") as fo:
            fo.write(undone)
        exit(1)
    else:
        print("Congratulation ! All your commonds done: {inputfile} ")
        logging.info( "Congratulation ! All your commonds done: {inputfile} ")
        os.system(f"touch {inputfile}.done") 



def check_done(file,handle="continue"):
    if os.path.exists(file):
        if handle == "break":
            print(f"# Attention cmd sikpped! due to {file} exisit")
            exit(0)
        else:
            logging.info(f"# Attention cmd sikpped! due to {file} exisit")
            return 1 


def cluster_run(cmds,cpus,queue,cluster="pbs"):

    pid = os.getpid()
    for cmd in cmds:
        if check_done("%s.done" %(cmd)) == 1 :
            pass
        else:
            cmd_basename = os.path.basename(cmd).replace(".sh","")
            commond=""
            if cluster == "pbs":
                commond='qsub -N %s_P%s -q %s -l nodes=1:ppn=%s -e %s.e -o %s.o %s' % (cmd_basename,pid,queue,cpus,cmd,cmd,cmd)
            if cluster == "slurm":
                commond='sbatch -J P%s_%s  -c %s  -e %s.e -o %s.o %s' %(pid,cmd_basename,cpus,cmd,cmd,cmd)
            else:
                commond='qsub -cwd -S /bin/bash -N P%s_%s -q %s -pe smp %s -e %s.e -o %s.o %s' % (pid,cmd_basename,queue,cpus,cmd,cmd,cmd)
            logging.info('Job submit: %s' % (commond))
            os.system(commond)

    keep(pid)


def keep(prefix):
    run_number=os.popen("squeue |grep %s|wc -l" % (prefix)).readlines()[0]
    while int(run_number) > 0 :
        #print int(run_number)
        time.sleep(20)
        run_number=os.popen("qstat |grep %s|wc -l" % (prefix)).read()
    return 1 



def local_run(cmd):
    stime = time.time()
    if check_done("%s.done" %(cmd)) == 1 :
        return 1 

    flag = os.system("sh %s  1>%s.std  2>%s.err && touch %s.done "% (cmd,cmd,cmd,cmd))
    if flag == 0 :
        logging.info("cmd Done!\t%s\t%s" % (spend_time(stime,time.time()),cmd))
    else :
        logging.error("cmd Fail!\t%s\t%s" % (spend_time(stime,time.time()),cmd))


def creat_sub_cmd(inputfile,line_num,name):

    job_list=[]
    os.system("mkdir -p %s.run" %(inputfile))
    task_list = open(inputfile,"r").read().splitlines()

    shell_num = int((len(task_list))/line_num)
    zfill_num = len(str(shell_num))

    for i in range(shell_num):
        num = i + 1 
        of=open("%s.run/%s%s.sh"%(inputfile,name, str(num).zfill(zfill_num) ),"w")
        of.write("#!/bin/bash\nset -e\ndate && hostname\n")
        of.write(f"cd {os.path.dirname(inputfile)}\n")
        for y in range(i*line_num,(i+1)*line_num):
            of.write(task_list[y]+"\n")
        of.write("date && touch %s.run/%s%s.sh.done\n" % (inputfile,name, str(num).zfill(zfill_num) ))
        of.close()
        job_list.append("%s.run/%s%s.sh"%(inputfile,name,str(num).zfill(zfill_num)))

    return job_list



def check_input(input):

    inputfile = os.path.abspath(input)

    if input == "std":
        with open("commond.sh","w") as of:
            for line in sys.stdin.readline():
                of.write(line)
        inputfile = os.getcwd()+"/commond.sh"

    return inputfile 


def read_cmd(input):
    '''
        read input form file or pipe and return a list of cmds
    '''
    if input == "std":
        return sys.stdin.read().splitlines()
    else:
        return open(input, 'r').read().splitlines()


def spend_time(stime,etime):
    spend=int(etime)-int(stime)
    h = int(spend/3600)
    remain = spend%3600
    m = int(remain/60)
    s = remain%60
    return "%s h %s m %s s" %(h,m,s)


def get_args():


    cpu_num=cpu_count()-2 
    parser = argparse.ArgumentParser(
    formatter_class = HelpFormatter,
    description = '''
RunCmd :
    Funciton:  Auto run  your mutiple task on local-host or SGE and PBS and Slurm cluster.
    Writer:    Meng Minghui < hyacinlee@163.com >
'''
    )
    parser.add_argument('-n',metavar='name',help='jobs cmd name  ',type=str,default="work")
    parser.add_argument('-l',metavar='line',help='how many line in one subprocess ',type=int,default=1)
    parser.add_argument('-f',help='froce to run the cmds without check done ',action="store_true")
    parser.add_argument('-t',metavar='type',help='set local, sge, slurm or pbs',default="local")
    parser.add_argument('-m',metavar='max',help='max number of local subprocess,defult is the number of cpus ',type=int,default=cpu_num)
    parser.add_argument('--ncpu',metavar='ncpu',help='number of cpu use for one task in sge | pbs ',default="2")
    parser.add_argument('--queue',metavar='queue',help='queue name for  sge | pbs ',default="all.q")
    parser.add_argument('INPUT',metavar='input',help='cmd input file must be given')

    args = parser.parse_args()

    return args



if __name__ == '__main__':
    main()
