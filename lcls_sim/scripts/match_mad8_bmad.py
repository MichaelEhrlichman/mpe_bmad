#!/bin/env python3

import sys

eps = 1e-6

bmad_file_name = sys.argv[1]
mad_file_name = sys.argv[2]

bmad_lines = []
with open(bmad_file_name,'r') as f:
  for line in f:
    if line.lstrip()[0] != '!':
      bmad_lines.append(line)

mad_lines = []
with open(mad_file_name,'r') as f:
  for line in f:
    if line.lstrip()[0] not in ['!','*','@','$']:
      mad_lines.append(line)

with open('mad.bmad_matched.twiss','w') as f:
  mad_ptr = 0
  for bmad_line in bmad_lines:
    bmad_line_data = bmad_line.split()
    s_bmad = float(bmad_line_data[1])
    while True:
      mad_line = mad_lines[mad_ptr]
      mad_line_data = mad_line.split()
      s_mad = float(mad_line_data[1])
      if abs(s_bmad-s_mad) < eps:
        f.write(mad_line)
        break
      mad_ptr += 1
