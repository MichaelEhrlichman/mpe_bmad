#!/bin/env python3

import sys

eps = 1e-6

bmad_file_name = sys.argv[1]
mad_file_name = sys.argv[2]

comments = ['!','*','@','$']

def parse_file(filename):
  data = []
  headers = []
  with open(filename,'r') as f:
    for line in f:
      if line.lstrip()[0] not in comments:
        data.append(line)
      else:
        headers.append(line)
  return data, headers

bmad_lines, bmad_headers = parse_file(bmad_file_name)
mad_lines, mad_headers = parse_file(mad_file_name)

with open('mad_bmad-matched.twiss','w') as f:
  mad_ptr = 0
  for header in mad_headers: 
    f.write(header)
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
