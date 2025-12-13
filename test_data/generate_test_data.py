#!/usr/bin/env python3
# Generate matrix.txt and act.txt for tests.
# matrix.txt: 1024 rows, each row 512 two-bit codes (00/01/10/11) separated by spaces
# act.txt: 512 complex numbers as int16_t, written as "re im" pairs (512 lines)

import os
import random

ROWS = 1024
COLS = 512
MATRIX_FILE = os.path.join(os.path.dirname(__file__), './', 'matrix.txt')
ACT_FILE = os.path.join(os.path.dirname(__file__), './', 'act.txt')

random.seed(42)
codes = ['00', '01', '10', '11']

# write matrix
with open(MATRIX_FILE, 'w', encoding='utf-8') as f:
    for r in range(ROWS):
        row_codes = [random.choice(codes) for _ in range(COLS)]
        f.write(' '.join(row_codes))
        f.write('\n')

# write act (512 complex numbers; each line `re im`)
with open(ACT_FILE, 'w', encoding='utf-8') as f:
    for i in range(COLS):
        re = 1
        im = 1
        f.write(f"{re} {im}\n")

print("Generated:", MATRIX_FILE, ACT_FILE)
# filepath: /Users/yipingwang/Documents/code/cpp/workspace/lut_c/scripts/generate_matrix_and_act.py
#!/usr/bin/env python3
# Generate matrix.txt and act.txt for tests.
# matrix.txt: 1024 rows, each row 512 two-bit codes (00/01/10/11) separated by spaces
# act.txt: 512 complex numbers as int16_t, written as "re im" pairs (512 lines)

import os
import random

ROWS = 1024
COLS = 512
MATRIX_FILE = os.path.join(os.path.dirname(__file__), '..', 'matrix.txt')
ACT_FILE = os.path.join(os.path.dirname(__file__), '..', 'act.txt')

random.seed(42)
codes = ['00', '01', '10', '11']

# write matrix
with open(MATRIX_FILE, 'w', encoding='utf-8') as f:
    for r in range(ROWS):
        row_codes = [random.choice(codes) for _ in range(COLS)]
        f.write(' '.join(row_codes))
        f.write('\n')

# write act (512 complex numbers; each line `re im`)
with open(ACT_FILE, 'w', encoding='utf-8') as f:
    for i in range(COLS):
        re = random.randint(-32768, 32767)
        im = random.randint(-32768, 32767)
        f.write(f"{re} {im}\n")

print("Generated:", MATRIX_FILE, ACT_FILE)