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
            value = random.choice([4.0, -4.0]) + random.uniform(-0.01, 0.01)
            if i % 2 == 0:
                value = random.choice([1, -1]) * 4000 + random.uniform(-10, 10)
            f.write(f"{value}\n")

if __name__ == '__main__':
    generate_act()
    generate_matrix()
