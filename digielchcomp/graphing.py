# -*- coding: utf-8 -*-
"""
Created on Mon Apr 24 10:11:17 2023

@author: Nathan
"""

import matplotlib.pyplot as plt
  
voltage = []
current = []
  
f = open('CV_kf_100_kb_2.0416.txt','r')
for row in f:
    row = row.split()
    if len(row) > 2:
        continue
    else:
        voltage.append(float(row[0]))
        current.append(float(row[1]))
  
plt.plot(voltage, current, color = 'g', label = 'File Data')
  
plt.cla()
plt.plot(voltage, current)
plt.title("'CV_kf_100_kb_2.0416.txt','r'")
plt.xlabel("Eapp [V]")
plt.ylabel("current [A]")