#!/usr/bin/env python3
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
import os

replacement_dict = {
    'cornua.uteri': 'cornua-uteri',
    'medulla.oblongata': 'medulla-oblongata',
    'optic.chiasm': 'optic-chiasm',
    'cerebral.cortex': 'cerebral-cortex',
    'corpus.uteri': 'corpus-uteri',
    'mammary.gland': 'mammary-gland',
    'bone.marrow': 'bone-marrow',
    'soft.horn': 'soft-horn',
    'lymph.node': 'lymph-node'
}

directory = '/vol2/zhangshiwen/sheep_cor/h3k27ac'
for filename in os.listdir(directory):
    if 'cpm' in filename and filename.endswith('.csv'):  #   'cpm'  CSVfile
        input_file = os.path.join(directory, filename)
        with open(input_file, 'r', newline='') as infile:
            lines = infile.readlines()

        if lines:  #  file
            modified_headers = []
            headers = lines[0].strip().split(',')
            for item in headers:
                for old, new in replacement_dict.items():
                    item = item.replace(old, new)
                modified_headers.append(item)
            modified_header = ','.join(modified_headers) + '\n'

        with open(input_file, 'w', newline='') as outfile:
            outfile.write(modified_header)  #
            outfile.writelines(lines[1:])  #  file

print("  'cpm'  file .")
