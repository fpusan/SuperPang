#!/usr/bin/env python3

MINIMAP2_PATH = 'minimap2'

import sys
from os.path import dirname, realpath
from os import mkdir
path = dirname(realpath(sys.argv[0]))
sys.path.remove(path)
sys.path.insert(0, realpath(path + '/../..'))
superpath = path + '/' + 'SuperPang.py'

import SuperPang
from inspect import getfile
data_path = getfile(SuperPang)
output_path = '/tmp/SuperPang_test_output'
datapath = dirname(realpath(data_path)) + '/' + 'test_data'

from subprocess import call
from glob import glob

print(f'\nAssuming that minimap2 can be found in the following path "{MINIMAP2_PATH}". Otherwise edit the MINIMAP2_PATH variable in this script\n')

ecode = call([superpath,'-f'] + glob(f'{datapath}/*fna') + ['--assume-complete', '-t', '12', '--force-overwrite', '-o', output_path, '--minimap2-path', MINIMAP2_PATH])

if not ecode:
    print(f'\n\nAll went well, I hope! Test results should be present in {output_path}\n')

