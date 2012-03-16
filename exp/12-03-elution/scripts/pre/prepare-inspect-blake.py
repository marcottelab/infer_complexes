#!/usr/bin/python
import os 
import sys
import stat
import mstb_helper as helper

usage_mesg = 'Usage: prepare-inspect.py'

dirname = 'inspect'

MSTB_HOME = helper.get_mstb_home()
CWD = helper.get_cwd()

conf = helper.read_conf( os.path.join(CWD, 'mstb.conf') )
if( len(conf) == 0 ):
    sys.exit(1)

helper.check_conf_file(conf,'DB_TRIE')
helper.check_conf_file(conf,'PATH_INSPECT')
helper.check_conf_file(conf,'DIR_INSPECT')

path_dir = os.path.join(CWD,dirname)
if( not os.access(path_dir,os.R_OK) ):
    sys.stderr.write('Create %s.\n'%(path_dir))
    os.mkdir(path_dir)

if( not os.path.isdir(path_dir) ):
    sys.stderr.write('\n%s is not a directory.\n'%path_dir)
    sys.stderr.write('Rename it and make %s directory.\n\n'%path_dir)
    sys.exit(1)

filename_cmd_tmpl = os.path.join(MSTB_HOME,'tmpl','inspect.cmd')
f_cmd_tmpl = open(filename_cmd_tmpl,'r')
cmd_tmpl = ''.join( f_cmd_tmpl.readlines() )
f_cmd_tmpl.close()

filename_in_tmpl = os.path.join(MSTB_HOME,'tmpl','inspect.in')
f_in_tmpl = open(filename_in_tmpl,'r')
in_tmpl = ''.join( f_in_tmpl.readlines() )
f_in_tmpl.close()

# Blake additions 3/11/2012
import numpy as np
n_per_run = 50
mzxmls_all = helper.get_mzxml_list()
nruns = np.floor(len(mzxmls_all) / n_per_run) + 1
for i in range(nruns):
    mzxmls = mzxmls_all[n_per_run*i:n_per_run*(i+1)]
    filename_run = os.path.join(CWD,'scripts','run-inspect-%s.sh' % i)
    f_run = open(filename_run,'w')
    f_run.write('#!/bin/bash\n')
    for basename_mzXML in mzxmls:
        filename_base = basename_mzXML.replace('.mzXML','')
        filename_in = os.path.join(CWD,dirname,'%s.inspect_in'%filename_base)
        filename_out = os.path.join(CWD,dirname,'%s.inspect_out'%filename_base)

        in_param = dict()
        in_param['DB_TRIE'] = conf['DB_TRIE']
        in_param['FILENAME_MZXML'] = os.path.join(CWD,'mzXML',basename_mzXML)

        cmd_param = dict()
        cmd_param['PATH_INSPECT'] = conf['PATH_INSPECT']
        cmd_param['DIR_INSPECT'] = conf['DIR_INSPECT']
        cmd_param['FILENAME_IN'] = filename_in
        cmd_param['FILENAME_OUT'] = filename_out
        f_run.write( cmd_tmpl.format(**cmd_param) )

        sys.stderr.write('Write %s.\n'%filename_in)
        f_in = open(filename_in,'w')
        f_in.write( in_tmpl.format(**in_param) )
        f_in.close()
    f_run.close()
    os.chmod(filename_run,stat.S_IRWXU)
sys.stderr.write('\nInsPecT is ready. Run %s.\n\n'%(filename_run))
