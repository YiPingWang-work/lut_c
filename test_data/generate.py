import os
import random

def generate_matrix():
    rows = 1024
    cols = 1024
    output_file = "./test_data/w.txt"
    with open(output_file, "w") as f:
        for _ in range(rows):
            row = [f"{random.randint(0, 1)}" for _ in range(cols*2)]
            # row = [f"00" for _ in range(cols)]
            # row = [f"01" for _ in range(cols)]
            # row = [f"10" for _ in range(cols)]
            # row = [f"11" for _ in range(cols)]
            f.write("".join(row) + "\n")

def generate_act():
    row = 2048
    output_file = "./test_data/act.txt"
    with open(output_file, "w") as f:
        for i in range(row):
            value = random.randint(-42, 42)
            # if i % 2 == 0:
            #     value /= 10
            # if value % 256 == 0:
            #     value = random.choice([-4.2, 4.2])
            # if value % 256 == 1:
            #     value = random.choice([-42, 42])
            f.write(f"{value}\n")

if __name__ == '__main__':
    generate_act()
    generate_matrix()
