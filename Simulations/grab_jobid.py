
import sys

with open(sys.argv[1],'r') as depend_file:
    lines = depend_file.readlines()
    lines = [line.rstrip() for line in lines]
    depend_id = lines[0].split('<')[1].split('>')[0]
    print(depend_id)

