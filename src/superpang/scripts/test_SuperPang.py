#!/usr/bin/env python3

MINIMAP2_PATH = 'minimap2'

import sys
from os.path import dirname, realpath
from os import mkdir
path = dirname(realpath(__file__))
if path in sys.path:
    sys.path.remove(path)
sys.path.insert(0, realpath(path + '/../..'))
superpath = path + '/' + 'SuperPang.py'

import superpang
from inspect import getfile
data_path = getfile(superpang)
output_path = '/tmp/SuperPang_test_output'
datapath = dirname(realpath(data_path)) + '/' + 'test_data'

from subprocess import call
from glob import glob

def main():
    print(f'\nAssuming that minimap2 can be found in the following path "{MINIMAP2_PATH}". Otherwise edit the MINIMAP2_PATH variable in this script\n')

    if len(sys.argv) > 1 and sys.argv[1] == 'fast':
        files = [f'{datapath}/2636415981.fna', f'{datapath}/2636416036.fna']
        iden = '-i1'
    elif len(sys.argv) > 1 and sys.argv[1] == 'fasth':
        files = [f'{datapath}/2636415981.fna', f'{datapath}/2636416036.fna']
        iden = '-i0.95'
    else:
       files = glob(f'{datapath}/*fna')
       iden = '-i0.95'

    args = [sys.executable, superpath,'-f'] + files + ['--assume-complete',
            iden, '-t12', '--force-overwrite', '-o', output_path, '--minimap2-path', MINIMAP2_PATH,
            '--keep-temporary', '--keep-intermediate', '--debug']
    argss = ' '.join(args)
    print(f'Running SuperPang with command "{argss}"')
    ecode = call(args)

    if not ecode:
        print(f'\n\nAll went well, I hope! Test results should be present in {output_path}\n')


if __name__ == '__main__':
    main()
