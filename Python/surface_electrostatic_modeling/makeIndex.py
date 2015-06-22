import os, sys

idxfile = "index.html"
iep_file = "93.iep"
all_files = os.listdir(os.getcwd())
ab_names = []
for fn in all_files:
    if fn[-4:] == '.pdb':
        ab_names.append(fn[:-4])

f1 = open(iep_file,"r")
iep = {}
idx =""
for line in f1.readlines():
    a = line.split()

    if len(a) == 7 and a[0]=='IEP': 
        idx = a[2]
        iep[idx] = []
    elif len(a) == 4 and a[0] == 'Isoelectric' and a[2] == '=': iep[idx].append(a[3])
    elif len(a) == 3 and a[0] == '7.00':
        iep[idx].append(a[2])
        idx = ''


media_root = 'bfx_analysis/20140806-1/'

code = """{% extends 'base_full.html' %}

          {% block title %}Electrostatic surface analysis results for LPAR1{% endblock %}

          {% block pagetitle %}<a href="index" class='brand'>Electrostatic surface analysis results for LPAR1</a>{% endblock %}
          {% load static %}
          {% block content %}
            <table class="table">
            <thead>
             <tr>
                 <th>XPA Number</th>
                 <th style="text-align: center;">Isoelectric point</th>
                 <th style="text-align: center;">Charge @ pH 7</th>
                 <th style="text-align: center;" width="250px">Left Side</th>
                 <th style="text-align: center;" width="250px">Right Side</th>
                 <th style="text-align: center;" width="250px">Top Side</th>
            </tr></thead>\n<tbody>"""

ab_names.sort()
for ab in ab_names:
    code += """<tr><th  width="100px">"""+ab+"""</th>
                   <td style="text-align: center;" width="100px">"""+iep[ab][0]+"""</td>
                   <td style="text-align: center;" width="100px">"""+iep[ab][1]+"""</td>
                   <td width="250px"><img src='{{ STATIC_URL }}"""+media_root+ab+"""_left.png'></td>
                   <td width="250px"><img src='{{ STATIC_URL }}"""+media_root+ab+"""_right.png'></td>
                   <td width="250px"><img src='{{ STATIC_URL }}"""+media_root+ab+"""_top.png'></td></tr>\n"""

code += "            </tbody></table> \n {% endblock %}"

f = open(idxfile,"w")
f.write(code)
f.close()