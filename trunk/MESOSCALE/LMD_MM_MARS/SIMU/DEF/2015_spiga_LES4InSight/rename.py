#! /usr/bin/env python
import os
path  = os.getcwd()
filenames = os.listdir(path)
for filename in filenames:
  print filename
  os.rename(filename, filename.replace(":", "-"))
