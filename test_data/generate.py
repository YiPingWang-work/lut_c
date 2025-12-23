import os
import random

def generate_matrix():
    rows = 10240
    cols = 10240
    output_file = "./test_data/w.txt"
    with open(output_file, "w") as f:
        for _ in range(rows):
            row = [f"{0}" for _ in range(cols*2)]
            f.write("".join(row) + "\n")

def generate_act():
    row = 20480
    output_file = "./test_data/act.txt"
    real = 0
    imag = 0
    with open(output_file, "w") as f:
        for i in range(row):
            value = random.uniform(-32, 31)
            if i%2 == 0:
                real += value
            else:
                imag += value
            # value = 1;
            f.write(f"{value}\n")
    print(f"Final accumulated complex value: {real} + {imag}i")

if __name__ == '__main__':
    generate_act()
    generate_matrix()
