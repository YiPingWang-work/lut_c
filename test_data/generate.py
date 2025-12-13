import os
import random

def generate_matrix():
    rows = 1024
    cols = 1024
    output_file = "./test_data/matrix_1024x1024.txt"
    with open(output_file, "w") as f:
        for _ in range(rows):
            row = [str(random.randint(0, 1)) for _ in range(cols)]
            f.write("".join(row) + "\n")

def generate_act():
    row = 2048
    output_file = "./test_data/act_1024.txt"
    with open(output_file, "w") as f:
        for _ in range(row):
            value = random.randint(-128, 127)
            f.write(f"{value}\n")

if __name__ == '__main__':
    generate_act()
